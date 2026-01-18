FROM ubuntu:22.04

# Avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    xz-utils \
    openbabel \
    nodejs \
    npm \
    && rm -rf /var/lib/apt/lists/*

# Download and install xTB
RUN wget -q https://github.com/grimme-lab/xtb/releases/download/v6.6.1/xtb-6.6.1-linux-x86_64.tar.xz \
    && tar -xf xtb-6.6.1-linux-x86_64.tar.xz \
    && mv xtb-6.6.1 /opt/xtb \
    && rm xtb-6.6.1-linux-x86_64.tar.xz

# Add xTB to PATH
ENV PATH="/opt/xtb/bin:${PATH}"
ENV XTBPATH="/opt/xtb/share/xtb"

# Create app directory
WORKDIR /app

# Copy package files
COPY package*.json ./

# Install Node.js dependencies
RUN npm install --production

# Copy application code
COPY . .

# Expose port
EXPOSE 3000

# Start the server
CMD ["node", "server.js"]
