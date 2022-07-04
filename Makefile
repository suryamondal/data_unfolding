
# LIB = -L${ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint  -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lMinuit -L/usr/lib64 -lc -lm -ldl -lcrypt -lpthread 

unfold:
	g++ -g -pthread -m64 -Wno-deprecated -std=c++11 -I${ROOTSYS}/include -o anal_muon_unfold anal_muon_unfold.cc `root-config --cflags --libs` -I/home/surya/products/RooUnfold/src/ -L/home/surya/products/RooUnfold/ -lRooUnfold

clean:
	rm anal_muon_unfold

unfold_basic:
	g++ -g -pthread -m64 -Wno-deprecated -std=c++11 -I${ROOTSYS}/include -o anal_muon_unfold_basic anal_muon_unfold_basic.cc `root-config --cflags --libs` -I/home/surya/products/RooUnfold/src/ -L/home/surya/products/RooUnfold/ -lRooUnfold

clean_basic:
	rm anal_muon_unfold_basic

unfold_basic_half:
	g++ -g -pthread -m64 -Wno-deprecated -std=c++11 -I${ROOTSYS}/include -o anal_muon_unfold_basic_half anal_muon_unfold_basic_half.cc `root-config --cflags --libs` -I/home/surya/products/RooUnfold/src/ -L/home/surya/products/RooUnfold/ -lRooUnfold

clean_basic_half:
	rm anal_muon_unfold_basic_halft

unfold_basic_SNM:
	g++ -g -pthread -m64 -Wno-deprecated -std=c++11 -I${ROOTSYS}/include -o anal_muon_unfold_SNM anal_muon_unfold_SNM.cc `root-config --cflags --libs` # -I/home/surya/products/RooUnfold/src/ -L/home/surya/products/RooUnfold/ -lRooUnfold

clean_basic_SNM:
	rm anal_muon_unfold_SNM
