# Get base image
FROM mathworks/matlab:r2022b

RUN ulimit -c unlimited

# Copy 2decomp files
COPY etc/2decomp_fft /usr/local/2decomp_fft

# Install dependencies
RUN sudo apt -y update
RUN sudo apt -y install make gfortran-9-multilib openmpi-bin openmpi-common openmpi-doc libopenmpi-dev

# Compile 2decomp
WORKDIR /usr/local/2decomp_fft
RUN sudo make all

# Go back to the main directory
WORKDIR /home/dns
ENTRYPOINT ["matlab"]