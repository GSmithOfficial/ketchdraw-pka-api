# Ketchdraw pKa API

FastAPI backend for pKa predictions, powering [Ketchdraw](https://gsmithofficial.github.io/Ketchdraw/).

## Features

- **Macro pKa**: Acidic and basic pKa values
- **Micro pKa**: Individual microstate equilibria  
- **logD**: pH-dependent lipophilicity
- **Microstate Distribution**: Population of each protonation state

## Powered By

- [UnipKa](https://github.com/finlayiainmaclean/unipka) by Finlay MacLean (Apache-2.0)
- Based on [Uni-pKa](https://github.com/dptech-corp/Uni-pKa) by DP Technology

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | API info |
| `/health` | GET | Health check |
| `/pka` | POST | Calculate pKa properties |
| `/docs` | GET | Interactive API documentation |

## Example Request
```bash
curl -X POST "https://your-app.railway.app/pka" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "pH": 7.4}'
```

## Local Development
```bash
pip install -r requirements.txt
uvicorn main:app --reload
```

## Deployment

Deployed on [Railway](https://railway.app).

## License

MIT License. UnipKa dependency is Apache-2.0.
