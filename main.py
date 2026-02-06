"""
Ketchdraw pKa API - FastAPI backend for pKa predictions.
Powered by UnipKa (Apache-2.0) by finlayiainmaclean.
"""
import asyncio
import traceback
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, field_validator
from typing import Optional

# Module-level import instead of per-request (Step 2.2)
from pka_logic import calculate_pka_properties, get_calculator


# ── Request model with input validation (Step 4.1) ──────────────────
class PKaRequest(BaseModel):
    smiles: str
    pH: Optional[float] = 7.4

    @field_validator("smiles")
    @classmethod
    def validate_smiles(cls, v):
        """Reject empty or excessively long SMILES before hitting the model."""
        v = v.strip()
        if not v:
            raise ValueError("SMILES string cannot be empty")
        # Guard against enormous molecules (proteins etc.) that would hang
        # 500 chars covers virtually all drug-like molecules
        if len(v) > 500:
            raise ValueError(
                f"SMILES too long ({len(v)} chars). "
                "Max 500 characters for drug-like molecules."
            )
        return v

    @field_validator("pH")
    @classmethod
    def validate_ph(cls, v):
        """Ensure pH is within a chemically meaningful range."""
        if v is not None and (v < 0 or v > 14):
            raise ValueError("pH must be between 0 and 14")
        return v


# ── Pre-warm the model on startup (Step 2.3) ────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    """Load the UnipKa model and run a warmup inference before accepting traffic."""
    print("⏳ Pre-warming UnipKa model...")
    await asyncio.to_thread(get_calculator)
    print("✅ UnipKa model loaded")
    
    # Run a throwaway prediction so PyTorch JIT-compiles everything
    # before real traffic arrives. Methane is the simplest molecule.
    print("⏳ Running warmup inference...")
    await asyncio.to_thread(calculate_pka_properties, "C", 7.4)
    print("✅ Warmup complete — ready for traffic")
    yield

# ── Create FastAPI app ───────────────────────────────────────────────
app = FastAPI(
    title="Ketchdraw pKa API",
    description="pKa prediction powered by UnipKa (Apache-2.0). "
                "Calculates macro/micro pKa, logD, and microstate distributions.",
    version="0.2.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan
)


# ── Global exception handler (Step 3.3 - no traceback in response) ──
@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    # Log full traceback to Railway logs for debugging
    print(f"UNHANDLED ERROR: {traceback.format_exc()}")
    # Return clean error to the client (no traceback exposure)
    return JSONResponse(
        status_code=500,
        content={
            "success": False,
            "error": str(exc)
        }
    )


# ── CORS ─────────────────────────────────────────────────────────────
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── Endpoints ────────────────────────────────────────────────────────
@app.get("/")
async def root():
    """API information and available endpoints."""
    return {
        "service": "Ketchdraw pKa API",
        "version": "0.2.0",
        "documentation": "/docs",
        "endpoints": {
            "/pka": "POST - Calculate pKa properties for a molecule",
            "/health": "GET - Health check"
        },
        "powered_by": "UnipKa by finlayiainmaclean (Apache-2.0)"
    }


@app.get("/health")
async def health():
    """
    Health check endpoint for Railway monitoring.
    Step 4.2: Verifies the model is actually loaded, not just that
    the web server is running.
    """
    model_loaded = get_calculator() is not None
    return {
        "status": "healthy" if model_loaded else "degraded",
        "model_loaded": model_loaded
    }


@app.post("/pka")
async def calculate_pka(request: PKaRequest):
    """
    Calculate pKa and related properties for a molecule.

    - **smiles**: SMILES string of the molecule (max 500 chars)
    - **pH**: pH value for distribution calculations (default: 7.4, range 0-14)
    """
    try:
        result = await asyncio.wait_for(
            asyncio.to_thread(calculate_pka_properties, request.smiles, request.pH),
            timeout=180
        )
        return result
    except asyncio.TimeoutError:
        return JSONResponse(
            status_code=504,
            content={
                "success": False,
                "error": "Prediction timed out (>120s). Molecule may be too complex.",
                "input_smiles": request.smiles
            }
        )
    except Exception as e:
        print(f"ERROR in /pka: {traceback.format_exc()}")
        return JSONResponse(
            status_code=500,
            content={
                "success": False,
                "error": str(e),
                "input_smiles": request.smiles
            }
        )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
