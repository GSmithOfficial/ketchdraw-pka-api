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
    
    Args:
        smiles: SMILES string of the molecule
        pH: pH value for calculations (default 7.4)
    
    Returns:
        Dictionary with pKa results
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
            "microstates": []
        }
    
    try:
        calc = get_calculator()
        
        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol)
        
        # Calculate pKa values
        acidic_pka = calc.get_acidic_macro_pka(mol)
        basic_pka = calc.get_basic_macro_pka(mol)
        
        # Get logD - NO pH parameter (uses internal default)
        try:
            logd = calc.get_logd(mol)
            logd = float(logd) if logd is not None and not math.isnan(logd) else None
        except Exception:
            logd = None
        
        # Get state penalty
        try:
            penalty_result = calc.get_state_penalty(mol, pH)
            if isinstance(penalty_result, tuple):
                state_penalty = float(penalty_result[0]) if not math.isnan(penalty_result[0]) else None
            else:
                state_penalty = float(penalty_result) if not math.isnan(penalty_result) else None
        except Exception:
            state_penalty = None
        
        # Get dominant microstate
        try:
            dominant = calc.get_dominant_microstate(mol, pH)
            dominant_smiles = Chem.MolToSmiles(dominant) if dominant else None
        except Exception:
            dominant_smiles = None
        
        # Get distribution
        microstates = []
        try:
            distribution_df = calc.get_distribution(mol, pH)
            if distribution_df is not None:
                for _, row in distribution_df.iterrows():
                    pop = row.get("population", 0)
                    microstates.append({
                        "smiles": str(row.get("smiles", "")),
                        "population": float(pop) if pop is not None and not math.isnan(pop) else 0
                    })
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
            "microstates": microstates
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
            "microstates": []
        }
