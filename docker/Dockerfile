FROM rust:1-bookworm AS builder

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        clang \
        cmake \
        git \
        libbz2-dev \
        libclang-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cargo build --release \
    --manifest-path rust/Cargo.toml \
    -p bioscript-cli

RUN install -m 0755 rust/target/release/bioscript /usr/local/bin/bioscript

FROM debian:bookworm-slim AS runtime

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libssl3 \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/bioscript /usr/local/bin/bioscript
RUN ln -s /usr/local/bin/bioscript /usr/local/bin/bs

LABEL org.opencontainers.image.source="https://github.com/OpenMined/bioscript"
LABEL org.opencontainers.image.description="BioScript runtime"
LABEL org.opencontainers.image.licenses="Apache-2.0"

WORKDIR /workspace

CMD ["bioscript"]
