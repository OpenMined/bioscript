FROM python:3.13-slim

LABEL org.opencontainers.image.source=https://github.com/OpenMined/bioscript
LABEL org.opencontainers.image.description="BioScript - Genetic variant classification tool"
LABEL org.opencontainers.image.licenses=Apache-2.0

# Install uv for fast package management
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set working directory
WORKDIR /app

# Copy Python package files
COPY python/pyproject.toml python/README.md ./
COPY python/src ./src

# Install bioscript
RUN uv pip install --system --no-cache .

# Set the default command to show help
CMD ["bioscript", "--help"]
