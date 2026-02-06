"""
pKa calculation logic using UnipKa
Lazy loading to avoid startup timeout
"""

import math

# Global calculator - loaded lazily
_calculator = None

def get_calculator():
    """Lazy load the UnipKa calculator"""
    global _calculator
    if _calculator is None:
        from unipka import UnipKa
        _calculator = UnipKa()
    return _calculator

def calculate_pka_properties(smiles: str, pH: float = 7.4) -> dict:
    """
    Calculate pKa properties for a molecule.
    """
    from rdkit import Chem
    
    # Validate SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": "Invalid SMILES string",
            "input_smiles": smiles,
            "canonical_smiles": None,
            "query_pH": pH,
            "pka": {},
            "logd": None,
            "state_penalty": None,
            "dominant_microstate": None,
            "microstates": [],
            "distribution_curve": None
        }
    
    try:
        calc = get_calculator()
        
        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol)
        
        # Calculate pKa values
        acidic_pka = calc.get_acidic_macro_pka(mol)
        basic_pka = calc.get_basic_macro_pka(mol)
        
        # Get logD - pH as KEYWORD argument
        try:
            logd = calc.get_logd(mol, pH=pH)
            logd = float(logd) if logd is not None and not math.isnan(logd) else None
        except Exception:
            logd = None
        
        # Get state penalty
        try:
            penalty_result = calc.get_state_penalty(mol, pH=pH)
            if isinstance(penalty_result, tuple):
                state_penalty = float(penalty_result[0]) if not math.isnan(penalty_result[0]) else None
            else:
                state_penalty = float(penalty_result) if not math.isnan(penalty_result) else None
        except Exception:
            state_penalty = None
        
        # Get dominant microstate
        try:
            dominant = calc.get_dominant_microstate(mol, pH=pH)
            dominant_smiles = Chem.MolToSmiles(dominant) if dominant else None
        except Exception:
            dominant_smiles = None
        
        # Get microstates at query pH
        microstates = []
        try:
            distribution_df = calc.get_distribution(mol, pH=pH)
            if distribution_df is not None:
                for _, row in distribution_df.iterrows():
                    pop = row.get("population", 0)
                    microstates.append({
                        "smiles": str(row.get("smiles", "")),
                        "charge": int(row.get("charge", 0)),
                        "population": float(pop) if pop is not None and not math.isnan(pop) else 0
                    })
        except Exception:
            pass
        
        # Get FULL distribution curve for plotting (pH 0-14)
        distribution_curve = None
        try:
            distribution_curve = get_distribution_curve(mol, calc)
        except Exception:
            pass
        
        return {
            "success": True,
            "input_smiles": smiles,
            "canonical_smiles": canonical_smiles,
            "query_pH": pH,
            "pka": {
                "acidic": float(acidic_pka) if acidic_pka is not None and not math.isnan(acidic_pka) else None,
                "basic": float(basic_pka) if basic_pka is not None and not math.isnan(basic_pka) else None
            },
            "logd": logd,
            "state_penalty": state_penalty,
            "dominant_microstate": dominant_smiles,
            "microstates": microstates,
            "distribution_curve": distribution_curve
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "input_smiles": smiles,
            "canonical_smiles": None,
            "query_pH": pH,
            "pka": {},
            "logd": None,
            "state_penalty": None,
            "dominant_microstate": None,
            "microstates": [],
            "distribution_curve": None
        }


def get_distribution_curve(mol, calc, ph_min=0, ph_max=14, steps=15):
    """
    Calculate distribution curve across pH range for plotting.
    Returns data structure ready for Chart.js or Plotly.
    
    OPTIMIZED: Single pass instead of double pass.
    Reduced default steps from 29 → 15 (every 1.0 pH unit) for 
    ~2x speedup with minimal visual difference in the curve.
    """
    import numpy as np
    from rdkit import Chem
    
    ph_values = np.linspace(ph_min, ph_max, steps).tolist()
    
    # Single pass: collect microstates AND populations together
    microstate_data = {}
    
    for i, ph in enumerate(ph_values):
        try:
            dist = calc.get_distribution(mol, pH=ph)
            if dist is not None:
                for _, row in dist.iterrows():
                    smi = str(row['smiles'])
                    
                    # First time seeing this microstate? Initialize it
                    if smi not in microstate_data:
                        microstate_data[smi] = {
                            "smiles": smi,
                            "charge": int(row['charge']),
                            "populations": [0.0] * len(ph_values)
                        }
                    
                    # Record population at this pH index
                    pop = row['population']
                    microstate_data[smi]["populations"][i] = (
                        float(pop) if not math.isnan(pop) else 0.0
                    )
        except Exception:
            pass
    
    return {
        "ph_values": ph_values,
        "microstates": list(microstate_data.values())
    }
