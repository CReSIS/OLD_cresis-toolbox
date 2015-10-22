stereo : stereo.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -DSTANDALONE -march=core2 -O3 -o stereo stereo.cpp -I . -lpng -ftree-vectorizer-verbose=1 

stereo_archive : stereo_archive.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -O3 -o stereo_archive stereo_archive.cpp -I . -lpng

stereo_archive2 : stereo_archive2.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -O3 -o stereo_archive2 stereo_archive2.cpp -I . -lpng

stereo_archive3 : stereo_archive3.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -O3 -o stereo_archive3 stereo_archive3.cpp -I . -lpng

compare2labeling : compare2labeling.cpp SImage.h SImageIO.h DTwoDimArray.h
	g++ -O3 -o compare2labeling compare2labeling.cpp -I . -lpng
