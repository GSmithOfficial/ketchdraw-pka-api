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
            "input_smiles": smiles
        }
    
    try:
        calc = get_calculator()
        
        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol)
        
        # Calculate pKa values
        acidic_pka = calc.get_acidic_macro_pka(mol)
        basic_pka = calc.get_basic_macro_pka(mol)
        
        # Get logD at specified pH
        logd = calc.get_logd(mol, pH)
        
        # Get state penalty
        penalty_result = calc.get_state_penalty(mol, pH)
        state_penalty = penalty_result[0] if isinstance(penalty_result, tuple) else penalty_result
        
        # Get dominant microstate
        dominant = calc.get_dominant_microstate(mol, pH)
        dominant_smiles = Chem.MolToSmiles(dominant) if dominant else None
        
        # Get distribution
        distribution_df = calc.get_distribution(mol, pH)
        microstates = []
        if distribution_df is not None:
            for _, row in distribution_df.iterrows():
                microstates.append({
                    "smiles": row.get("smiles", ""),
                    "population": float(row.get("population", 0)) if not math.isnan(row.get("population", 0)) else 0
                })
        
        return {
            "success": True,
            "input_smiles": smiles,
            "canonical_smiles": canonical_smiles,
            "query_pH": pH,
            "pka": {
                "acidic": float(acidic_pka) if not math.isnan(acidic_pka) else None,
                "basic": float(basic_pka) if not math.isnan(basic_pka) else None
            },
            "logd": float(logd) if not math.isnan(logd) else None,
            "state_penalty": float(state_penalty) if not math.isnan(state_penalty) else None,
            "dominant_microstate": dominant_smiles,
            "microstates": microstates
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "input_smiles": smiles
        }
