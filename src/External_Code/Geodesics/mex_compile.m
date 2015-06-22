mex ./geodesics/comp_geodesics_to_all.cpp
mex ./geodesics/comp_geodesics_pairs.cpp


% Added the clause '-Dchar16_t=uint16_t' to work with Clang coming with OS 5+.
% mex -Dchar16_t=uint16_t ./geodesics/comp_geodesics_to_all.cpp
% mex -Dchar16_t=uint16_t ./geodesics/comp_geodesics_pairs.cpp

% clr
% mex -DCHAR16_T=uint16_t ./geodesics/comp_geodesics_to_all.cpp
% mex -DCHAR16_T=uint16_t ./geodesics/comp_geodesics_pairs.cpp

% mex -Dchar16_t=UINT16_T ./geodesics/comp_geodesics_to_all.cpp
% mex -Dchar16_t=UINT16_T ./geodesics/comp_geodesics_pairs.cpp