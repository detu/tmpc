## Installing dependencies
- *CMake 3.9.0 or higher*. Note that Ubuntu 17.04 comes with CMake 3.7, therefore you may need to install latest CMake manually from https://cmake.org/download/
- *Boost* `sudo apt install libboost-dev`.
- *Eigen3 3.3.3 or higher*. Ubuntu 17.04 comes with Eigen3 3.3.2-1, therefore an up-to-date Eigen3 must be downloaded from http://eigen.tuxfamily.org and installed so that CMake can find it.
- *Blaze 3.3 or higher* https://bitbucket.org/blaze-lib/blaze.
- *BLASFEO* https://github.com/giaf/blasfeo (optional, only if one of `TMPC_WITH_HPMPC`, `TMPC_WITH_HPIPM` is selected). Select a proper target architecture by setting the `TARGET` variable in `Makefile.rule` or in `CMake`. Build and install as usual. The build system searches for BLASFEO in `/opt/blasfeo` by default.
- *HPMPC* https://github.com/giaf/hpmpc (optional, only if `TMPC_WITH_HPMPC` is selected). Build and install from sources. Make sure that in `Makefile.rule` `USE_BLASFEO` is set to `1` and `BLASFEO_PATH` is set to the correct path to BLASFEO installed on previous step. I experienced HPMPC crashing or giving incorrect results if built without BLASFEO. The build system searches for HMPMC in `/opt/hpmpc` by default.
- *HPIPM* https://github.com/giaf/hpipm (optional, only if `TMPC_WITH_HPMPC` is selected). Build and install from sources. Make sure that in `Makefile.rule` `BLASFEO_PATH` is set to the correct path to BLASFEO. The build system searches for HMIPM in `/opt/hpipm` by default.
- *qpOASES* https://projects.coin-or.org/qpOASES (optional, only if `TMPC_WITH_qpOASES` is selected). Build and install as usual.
- *CasADi 3.2.0 or higher* https://github.com/casadi/casadi/wiki/InstallationInstructions (optional, only if `TMPC_WITH_CASADI` is selected). Python bindings for CasADi must be installed. It is recommended to use Python3 CasADi bindings; I didn't test the latest code with Python2. To build Python3 CasADi bindings, specify `-DWITH_PYTHON3=ON` in `cmake` command line when configuring CasADi.
- *Google Test* https://github.com/google/googletest must be installed and findable by the CMake build system (optional, only if `TMPC_WITH_TEST` is selected).
- *Google Benchmark* https://github.com/google/benchmark must be installed and findable by the CMake build system (optional, only if `TMPC_WITH_BENCHMARK` is selected).
- *JSON for Modern C++*. Build and install from sources https://github.com/nlohmann/json (optional, only if `TMPC_WITH_JSON` is selected).

## Building
1. Install the dependencies.
2. Assuming that you are in the `tmpc` source root, do

    ```bash
    mkdir build && cd build
    ```
3. Run CMake
```bash
cmake -DTMPC_WITH_qpOASES=ON -DTMPC_WITH_HPMPC=ON -DTMPC_WITH_HPIPM=ON ..
```
4. Build
```bash
make -j 10
```
5. Run tests
```bash
ctest
```