"""
Ketchdraw pKa API - FastAPI backend for pKa predictions.
Powered by UnipKa (Apache-2.0) by finlayiainmaclean.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional

from pka_logic import calculate_pka_properties


# Request/Response models
class PKaRequest(BaseModel):
    smiles: str
    pH: Optional[float] = 7.4


class PKaResponse(BaseModel):
    input_smiles: str
    canonical_smiles: Optional[str] = None
    query_pH: float
    success: bool
    error: Optional[str] = None
    pka: dict
    logd: Optional[float] = None
    state_penalty: Optional[float] = None
    dominant_microstate: Optional[str] = None
    microstates: list


# Create FastAPI app
app = FastAPI(
    title="Ketchdraw pKa API",
    description="pKa prediction powered by UnipKa (Apache-2.0). "
                "Calculates macro/micro pKa, logD, and microstate distributions.",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# Enable CORS for Ketchdraw frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://gsmithofficial.github.io",  # Your GitHub Pages
        "http://localhost:3000",              # Local development
        "http://localhost:8080",
        "*"  # Remove in production if you want to restrict
    ],
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


@app.post("/pka", response_model=PKaResponse)
async def calculate_pka(request: PKaRequest):
    """Calculate pKa properties for a molecule"""
    import traceback
    try:
        result = calculate_pka_properties(request.smiles, request.pH)
        return result
    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"ERROR in /pka: {error_msg}")  # This will show in Railway logs
        return {
            "success": False,
            "error": str(e),
            "traceback": error_msg,
            "input_smiles": request.smiles
        }
    """
    Calculate pKa and related properties for a molecule.
    
    - **smiles**: SMILES string of the molecule
    - **pH**: pH value for distribution calculations (default: 7.4)
    
    Returns:
    - Acidic and basic macro pKa values
    - logD at the specified pH
    - State penalty (for permeability estimation)
    - Dominant microstate SMILES at the specified pH
    - Full list of microstates with populations
    """
    result = calculate_pka_properties(request.smiles, request.pH)
    return result


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
