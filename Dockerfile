FROM python:3.11-slim

# Install system dependencies required by RDKit
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libsm6 \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first (for caching)
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Expose port (documentation only)
EXPOSE 8000

# Run with shell to expand $PORT variable
CMD uvicorn main:app --host 0.0.0.0 --port ${PORT:-8000} --workers 2
