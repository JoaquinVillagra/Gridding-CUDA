default : gridding.cu 
	nvcc -Xptxas="-v" gridding.cu -o gridding
