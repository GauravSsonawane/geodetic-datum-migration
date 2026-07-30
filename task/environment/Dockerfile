# Single image, built for both the agent run and the verifier. Ground truth is
# overlaid at /tests only at verify time, so it is never readable during the
# agent's run.
FROM public.ecr.aws/docker/library/ubuntu:24.04@sha256:0d39fcc8335d6d74d5502f6df2d30119ff4790ebbb60b364818d5112d9e3e932

# Runtime and verifier dependencies, baked at build time (pinned) so the verifier
# installs nothing at verify time.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip ca-certificates \
    && pip install --no-cache-dir --break-system-packages \
        numpy==2.3.5 pytest==8.4.1 pytest-json-ctrf==0.3.5 \
    && rm -rf /var/lib/apt/lists/*

# The legacy archive the agent has to migrate.
COPY data /app/data

# Parent directory of the declared artifacts.
RUN mkdir -p /app/out

WORKDIR /app
