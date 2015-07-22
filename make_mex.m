% Debug
mex -g -DPLATFORM_TYPE=rtmc::CyberMotion -DN_PREDICTION=20 -DS_FUNCTION_NAME=controller_cms -output controller_cms controller.cpp build\src\Debug\rtmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen

% Release
mex -DPLATFORM_TYPE=rtmc::CyberMotion -DN_PREDICTION=20 -DS_FUNCTION_NAME=controller_cms -output controller_cms controller.cpp build\src\Release\rtmc.lib E:\software\qpOASES\build\libs\Release\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen

mex -g -DPLATFORM_TYPE=CyberMotion -DS_FUNCTION_NAME=plant_output_cms -output plant_output_cms plant_output.cpp build\src\Debug\rtmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen
mex -g -DPLATFORM_TYPE=MotionPlatformX -DN_PREDICTION=2 -DS_FUNCTION_NAME=controller_x -output controller_x controller.cpp build\src\Debug\rtmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen
mex -g -DPLATFORM_TYPE=MotionPlatformX -DS_FUNCTION_NAME=plant_output_x -output plant_output_x plant_output.cpp build\src\Debug\rtmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\ 