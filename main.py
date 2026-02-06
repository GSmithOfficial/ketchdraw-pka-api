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
from pydantic import BaseModel
from typing import Optional

# Module-level import instead of per-request (Step 2.2)
from pka_logic import calculate_pka_properties, get_calculator


# Request model
class PKaRequest(BaseModel):
    smiles: str
    pH: Optional[float] = 7.4


# Pre-warm the model on startup so first user doesn't wait (Step 2.3)
@asynccontextmanager
async def lifespan(app: FastAPI):
    """Load the UnipKa model before accepting traffic."""
    print("⏳ Pre-warming UnipKa model...")
    await asyncio.to_thread(get_calculator)
    print("✅ UnipKa model ready")
    yield


# Create FastAPI app
app = FastAPI(
    title="Ketchdraw pKa API",
    description="pKa prediction powered by UnipKa (Apache-2.0). "
                "Calculates macro/micro pKa, logD, and microstate distributions.",
    version="0.2.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan
)

# Global exception handler - catches ALL errors and returns details
@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    error_msg = traceback.format_exc()
    print(f"UNHANDLED ERROR: {error_msg}")
    return JSONResponse(
        status_code=500,
        content={
            "success": False,
            "error": str(exc),
            "traceback": error_msg
        }
    )

# Enable CORS for Ketchdraw frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def root():
    """API information and available endpoints."""
    return {
        "service": "Ketchdraw pKa API",
        "version": "0.1.0",
        "documentation": "/docs",
        "endpoints": {
            "/pka": "POST - Calculate pKa properties for a molecule",
            "/health": "GET - Health check"
        },
        "powered_by": "UnipKa by finlayiainmaclean (Apache-2.0)"
    }


@app.get("/health")
async def health():
    """Health check endpoint for Railway monitoring."""
    return {"status": "healthy"}


import asyncio

@app.post("/pka")
async def calculate_pka(request: PKaRequest):
    """
    Calculate pKa and related properties for a molecule.
    
    - **smiles**: SMILES string of the molecule
    - **pH**: pH value for distribution calculations (default: 7.4)
    """
   try:
        # No more "from pka_logic import" here — already imported at top
        result = await asyncio.wait_for(
            asyncio.to_thread(calculate_pka_properties, request.smiles, request.pH),
            timeout=120
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
        error_msg = traceback.format_exc()
        print(f"ERROR in /pka: {error_msg}")
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
