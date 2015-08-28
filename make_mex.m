% Debug
mex -g -DPLATFORM_TYPE=mpmc::CyberMotion -DN_PREDICTION=20 -DS_FUNCTION_NAME=controller_cms -output controller_cms controller.cpp build\lib\Debug\mpmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen

% Release
mex -g -DPLATFORM_TYPE=mpmc::CyberMotion -DN_PREDICTION=20 -DS_FUNCTION_NAME=controller_cms -output controller_cms controller.cpp build\lib\Release\mpmc.lib E:\software\qpOASES\build\libs\Release\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen

mex -g -DPLATFORM_TYPE=mpmc::CyberMotion -DS_FUNCTION_NAME=plant_output_cms -output plant_output_cms plant_output.cpp build\lib\Debug\mpmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen

% X - debug
mex -g -DPLATFORM_TYPE=mpmc::MotionPlatformX -DN_PREDICTION=2 -DS_FUNCTION_NAME=controller_x -output controller_x controller.cpp build\lib\Debug\mpmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen
% X - release
mex -g -DPLATFORM_TYPE=mpmc::MotionPlatformX -DN_PREDICTION=2 -DS_FUNCTION_NAME=controller_x -output controller_x controller.cpp build\lib\Release\mpmc.lib E:\software\qpOASES\build\libs\Release\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\eigen
% mex -g -DPLATFORM_TYPE=MotionPlatformX -DS_FUNCTION_NAME=plant_output_x -output plant_output_x plant_output.cpp build\lib\Debug\mpmc.lib E:\software\qpOASES\build\libs\Debug\qpOASES.lib -DWIN32 -Iinclude -IE:\software\qpOASES\include -IE:\software\ 