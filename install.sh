sudo apt-get install -y libgl1-mesa-dev libglu1-mesa-dev libx11-dev gfortran libfontconfig-dev curl zip unzip tar
./external/vcpkg/bootstrap-vcpkg.sh
./external/vcpkg install cgal OpenCASCADE glog Ceres argparse Geogram