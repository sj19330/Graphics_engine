#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <ModelTriangle.h>
#include <fstream>
#include <sstream>
#include "glm/ext.hpp"
#include <unordered_map>
#include <RayTriangleIntersection.h>
#include <cmath>
#include <TextureMap.h>


#define WIDTH 320*2
#define HEIGHT 240*2
using namespace std;
bool Quit = false;

glm::vec3 CAMERAPOSITION = 	glm::vec3(0.0, 0.0, 4.0);
glm::mat3 CAMERAORIENTATION = glm::mat3(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);

glm::vec3 LIGHTPOSITION = glm::vec3(0.0, 0.4, 0.0);
vector<glm::vec3> LIGHTPOSITIONS;

void initiliseLightPositions(){
	LIGHTPOSITIONS.clear();
	LIGHTPOSITIONS.push_back(LIGHTPOSITION);
	for (float i=0.01; i<0.08; i=i+0.01){
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(i,0.0,0.0));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(-i,0.0,0.0));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(0.0,0.0,i));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(0.0,0.0,-i));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(i,0.0,i));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(i,0.0,-i));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(-i,0.0,i));
		LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(-i,0.0,-i));
		//  LIGHTPOSITIONS.push_back(LIGHTPOSITION + glm::vec3(0.0,i,0.0));
	}
}

//								************################# SECTION FOR PARSER #################************
struct getNumbersFromLineReturn{
	glm::vec3 coords;
	glm::vec3 normals;
	glm::vec3 textures;
}typedef getNumbersFromLineReturn;

getNumbersFromLineReturn getNumbersFromLine(string str){
	float coord;
	vector<float> coords;
	vector<float> textures;
	vector<float> normals;
	glm::vec3 empty;
	string emptyString;
	vector<string> seperatedLine = split(str, ' ');
    for (int i=0; i<seperatedLine.size(); i++){
		vector<string> seperatedWord = split(seperatedLine[i], '/');
		if((seperatedWord.size() == 1)){
			if(stringstream(seperatedLine[i]) >> coord){
				coords.push_back(coord);
			}
		}
		if(seperatedWord.size() == 2){
            coords.push_back(stoi(seperatedWord[0]));
            if(seperatedWord[1] != emptyString){
                textures.push_back(stoi(seperatedWord[1]));
            }
        }
		if(seperatedWord.size() == 3){
			coords.push_back(stoi(seperatedWord[0]));
			normals.push_back(stoi(seperatedWord[2]));
		}
    }
	getNumbersFromLineReturn returnValues;
	returnValues.coords = glm::vec3(coords[0],coords[1],coords[2]);
	if(textures.size() == 3){
        returnValues.textures = glm::vec3(textures[0], textures[1], textures[2]);
    }//may want an else here
	if (normals.size() == 3){
		returnValues.normals = glm::vec3(normals[0], normals[1], normals[2]);
	}else returnValues.normals = empty;
	
	
	return returnValues;
}


string getWordFromLine(string str){
	stringstream line;
	string word;

	line << str;
	while(!line.eof()){
		line >> word;
	}
	return word;

}
//read in both mtl and obj file and return a vector of model triangles with their colour
vector<ModelTriangle> parser(){
	float scale = 0.17;
	string line;
	string cline;
	vector<glm::vec3> vertexes;
	vector<glm::vec3> vertexNormals;
	vector<TexturePoint> vertexTextures;
	vector<glm::vec3> empty;
	glm::vec3 empty2;
	vector<glm::vec3> flines;
	vector<ModelTriangle> triangles;
	ifstream file;
	ifstream cfile;
	unordered_map<string, glm::vec3> colours;
	string colourForThisObject = " ";
	string isThisObjectMirrored = "n";
	bool mirrored;

	//open read and creat a hashmap of mtl file contents then close
	cfile.open("src/materials.mtl");
	if(!cfile.is_open()){
		cout << "error opening mtl file" << endl;
	}
	else{
		while(cfile.good()){
			getline(cfile, cline);
			if(cline[0]=='n'){
				string tempColour = getWordFromLine(cline);
				getline(cfile, cline);
				glm::vec3 rgb = getNumbersFromLine(cline).coords;
				colours[tempColour] = rgb;
			}
		}
	}
	cfile.close();
	//open and read the contents of the obj file then close
	file.open("src/sphere2.obj");
	if (!file.is_open()){
		cout << "error opening obj file" << endl;
	}
	else{
		while(file.good()){
			getline(file, line);
			if(line[0] == 'u'){
				colourForThisObject = getWordFromLine(line);
			}
			else if(line[0] == 'v' && line[1] == ' '){
				glm::vec3 vertex = getNumbersFromLine(line).coords;
				vertexes.push_back(vertex);
			}
			else if(line[0] == 'v' && line[1] == 'n'){
				glm::vec3 vertexNormal = getNumbersFromLine(line).coords;
				vertexNormals.push_back(vertexNormal);
			}
			else if(line[0] == 'v' && line[1] == 't'){
                glm::vec3 temp = getNumbersFromLine(line).coords;
                TexturePoint vertexTexture = TexturePoint(temp[0], temp[1]);
                vertexTextures.push_back(vertexTexture);
            }
			else if(line[0] == 'm'){
				isThisObjectMirrored = getWordFromLine(line);
			}
			else if(line[0] == 'f'){
				getNumbersFromLineReturn fline = getNumbersFromLine(line);
				Colour col = Colour(colours[colourForThisObject][0]*255, colours[colourForThisObject][1]*255, colours[colourForThisObject][2]*255);
				if (isThisObjectMirrored == "y"){
					mirrored = true;
				}
				else{
					mirrored = false;
				}
				ModelTriangle temp = ModelTriangle(vertexes[((fline.coords[0])-1)]*scale,vertexes[((fline.coords[1])-1)]*scale,vertexes[((fline.coords[2])-1)]*scale, col);
				temp.mirrored = mirrored;
				if (vertexNormals != empty){
					if(fline.normals == empty2){
						temp.vertexNormals = {vertexNormals[(fline.coords[0])-1], vertexNormals[(fline.coords[1])-1], vertexNormals[(fline.coords[2])-1]};	
					}
					else{
						temp.vertexNormals = {vertexNormals[(fline.normals[0])-1], vertexNormals[(fline.normals[1])-1], vertexNormals[(fline.normals[2])-1]};
					}
				}
				if(vertexTextures.size() != 0){
                    if(fline.textures == empty2){
                    }else{
                        temp.texturePoints = {vertexTextures[(fline.textures[0])-1], vertexTextures[(fline.textures[1])-1], vertexTextures[(fline.textures[2])-1]};
                        temp.textured = true;
                    }
                    
                }
				triangles.push_back(temp);
			}
		}
	}
	file.close();
	return triangles;
}



//								************################# SECTION FOR RASTURISING METHOD #################************
float DEPTHBUFFER[WIDTH][HEIGHT];

void initialiseDepthBuffer(){
	for(int i=0; i<WIDTH; i++){
		for(int j=0; j<HEIGHT; j++){
			DEPTHBUFFER[i][j] = 0;
		}
	}
}


CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength){
	float scaleUp = 600;
	glm::vec3 cameraToVertex = vertexPosition - CAMERAPOSITION;
	glm::vec3 adjusted = cameraToVertex*CAMERAORIENTATION;
	float depth = adjusted.z;
	float canvasX = (focalLength*((adjusted.x)/depth))*scaleUp +WIDTH/2;
	float canvasY = (focalLength*((adjusted.y)/depth))*scaleUp +HEIGHT/2;
	canvasX = WIDTH - canvasX;
	return CanvasPoint(canvasX, canvasY, depth);
}
void drawPoint(DrawingWindow &window, CanvasPoint point) {
	uint32_t white = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	window.setPixelColour(point.x, point.y, white);
}
void drawLineConsideringDepth(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour col){
	uint32_t colour = (255<<24) + (int(col.red) << 16) +(int(col.green) << 8) +int(col.blue);
	float xdiff = to.x - from.x;
	float ydiff = to.y - from.y;
	float depthDiff = to.depth - from.depth;
	float numberOfSteps = fmax(abs(xdiff), abs(ydiff));
	float individualStepDepthDifference = depthDiff/numberOfSteps;
	float xStepSize = xdiff/numberOfSteps;
	float yStepSize = ydiff/numberOfSteps;
	for (float i=0.0; i<numberOfSteps; i++){
		int x = round(from.x+(xStepSize*i));
		int y = from.y+(yStepSize*i);
		float depth = from.depth+individualStepDepthDifference*i;
		if(DEPTHBUFFER[x][y] == 0){
			window.setPixelColour(round(x), round(y), colour);
			DEPTHBUFFER[x][y] = 1/depth;
		}
		else if(DEPTHBUFFER[x][y]>(1/depth)){
			window.setPixelColour(round(x), round(y), colour);
			DEPTHBUFFER[x][y] = 1/depth;
		}
		
	}
}
void drawEmptyTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour col){
	drawLineConsideringDepth(window, triangle.v1(), triangle.v0(), col);
	drawLineConsideringDepth(window, triangle.v1(), triangle.v2(), col);
	drawLineConsideringDepth(window, triangle.v2(), triangle.v0(), col);
}
CanvasTriangle sort(CanvasTriangle tri){
	CanvasPoint top = tri.v0();
	CanvasPoint middle = tri.v1();
	CanvasPoint bottom = tri.v2();
	if(top.y < middle.y){
		CanvasPoint temp = top;
		top = middle;
		middle = temp;
	}
	if(top.y < bottom.y){
		CanvasPoint temp = top;
		top = bottom;
		bottom = temp;
	}
	if(middle.y < bottom.y){
		CanvasPoint temp = middle;
		middle = bottom;
		bottom = temp;
	}
	CanvasTriangle sorted = CanvasTriangle(top ,middle, bottom);
	return sorted;
}
CanvasPoint findNewPoint(CanvasTriangle tri){
	float heightDiff = tri.v0().y - tri.v2().y;
	float widthDiff = tri.v0().x - tri.v2().x;
	float depthDiff = tri.v0().depth - tri.v2().depth;
	float ratio1 = -(widthDiff/heightDiff);
	//line below could be -ve, not too sure
	float ratio2 = -(depthDiff/heightDiff);
	float heighDiff2 = tri.v0().y - tri.v1().y;
	float newz = tri.v0().depth + ratio2*heighDiff2;
	float newx = tri.v0().x + ratio1*heighDiff2;
	float newy = tri.v1().y;


	CanvasPoint result = CanvasPoint(newx, newy, newz);
	return result;
}
vector<CanvasPoint> interpolateLine(CanvasPoint start, CanvasPoint end, int side){
	vector<CanvasPoint> v;
	int numberOfValues = abs(start.y - end.y);
	float individualSeperation = (start.x - end.x)/numberOfValues;
	//may want to check sign on line below
	float individualDepthSepperation = -(start.depth - end.depth)/numberOfValues;
	if(side == 0){
		for(int i=0; i<=numberOfValues; i++){ 
			v.push_back(CanvasPoint(roundf(start.x + i*(-individualSeperation)), start.y - i, (start.depth + i*(individualDepthSepperation))));
		}
	}
	else if(side == 1){
		for(int i=0; i<=numberOfValues; i++){ 
			v.push_back(CanvasPoint(roundf(start.x + i*(-individualSeperation)), start.y + i, (start.depth + i*(individualDepthSepperation))));
		}
	}
	return v;
}
void interpolateSidesAndFill(CanvasTriangle tri, DrawingWindow &window, Colour col, int side){
	vector<CanvasPoint> side1 = interpolateLine(tri.v0(), tri.v1(), side);
	vector<CanvasPoint> side2 = interpolateLine(tri.v0(), tri.v2(), side);
	for(int i = 0; i< side1.size(); i++){
		drawLineConsideringDepth(window, side1[i], side2[i], col);
	}
	drawEmptyTriangle(window, tri, col);
}
void drawFilledTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col){
	//sort the points in terms of height 
	tri = sort(tri);
	//find the coord of middle point on the side
	CanvasPoint newp = findNewPoint(tri);
	//make two triangles out of the main one 
	CanvasTriangle topTri = CanvasTriangle(tri.v0(), tri.v1(), newp);
	CanvasTriangle bottomTri = CanvasTriangle(tri.v2(), tri.v1(), newp);
	//fill each triangle with the rasturising method
		//interpolate each side
	interpolateSidesAndFill(topTri, window, col, 0);
	interpolateSidesAndFill(bottomTri, window, col, 1);
}


//     texture map section



vector<CanvasPoint> interpolateSideForTex(CanvasPoint start, CanvasPoint end){
    vector<CanvasPoint> v;
	// cout << start.y << " " << end.y << endl;
    float numberOfValuesy = abs(round(end.y) - round(start.y));
	// cout << numberOfValuesy <<endl;
    float numberOfValuesx = abs(round(end.x) - round(start.x));
    float numberOfValues = fmax(numberOfValuesy, numberOfValuesx);
    float stepSizex = (end.x - start.x)/numberOfValues;
    float stepSizey = (end.y - start.y)/numberOfValues;
    float individualDepthSepperation = -(start.depth - end.depth)/numberOfValues;
    for (int i=0; i<numberOfValues; i++){
        CanvasPoint temp = CanvasPoint(round(start.x + i*stepSizex), round(start.y + i*stepSizey), (start.depth + i*(individualDepthSepperation)));
        v.push_back(temp);
    }
    return v;
}

vector<CanvasPoint> interpolateSideForTex(TexturePoint start, TexturePoint end){
    vector<CanvasPoint> v;
    float numberOfTexValuesy = abs(round(end.y) - round(start.y));
	// cout << numberOfTexValuesy << endl;
    float numberOfTexValuesx = abs(round(end.x) - round(start.x));
	// cout << numberOfTexValuesx << endl;
    float numberOfTexValues = fmax(numberOfTexValuesy, numberOfTexValuesx);
    float texStepSizex = (end.x - start.x)/numberOfTexValues;
    float texStepSizey = (end.y - start.y)/numberOfTexValues;
    for (int i=0; i<numberOfTexValues; i++){
        CanvasPoint temp = CanvasPoint(round(start.x + i*texStepSizex), round(start.y + i *texStepSizey));
        v.push_back(temp);
    }
    return v;
}

CanvasPoint findNewPointTex(CanvasTriangle tri){
	float heightDiff = tri.v0().y - tri.v2().y;
	float widthDiff = tri.v0().x - tri.v2().x;
	float depthDiff = tri.v0().depth - tri.v2().depth;
	float proportionDownSide = (tri.v0().y - tri.v1().y)/heightDiff;
	float newz = tri.v0().depth - depthDiff*proportionDownSide;
    float newx = tri.v0().x - widthDiff*proportionDownSide;
    float newy = tri.v1().y;
    
    float texHeightDiff = tri.v0().texturePoint.y - tri.v2().texturePoint.y;
    float texWidthDiff = tri.v0().texturePoint.x - tri.v2().texturePoint.x;
    float newTexX = tri.v0().texturePoint.x - texWidthDiff*proportionDownSide;
    float newTexY = tri.v0().texturePoint.y - texHeightDiff*proportionDownSide;
	
    TexturePoint newTex = TexturePoint(round(newTexX), round(newTexY));
	CanvasPoint result = CanvasPoint(round(newx), newy, newz);
    result.texturePoint = newTex;
	return result;
}


void drawTextureLine(CanvasPoint texStart, CanvasPoint texEnd, CanvasPoint from, CanvasPoint to, TextureMap map, DrawingWindow &window){
    // cout << " in drawTextureline" << endl;
	// cout << map.width <<" x " << map.height << endl;
	// cout << from.x << " - " << to.x <<endl;
	float numberOfSteps = abs(to.x - from.x);
	// cout << texStart << " " << texEnd << endl;
    float depthDiff = to.depth - from.depth;
	float individualStepDepthDifference = depthDiff/numberOfSteps;
    vector<CanvasPoint> textureLine = interpolateSideForTex(texStart, texEnd);
	// cout << textureLine.size() <<  " " << numberOfSteps << endl;
    for(float i=0.0; i<numberOfSteps; i++){
        float oneStepThroughTexLine = (textureLine.size()-1)/numberOfSteps;
        int x = from.x - i;
        int y = from.y;
        float depth = from.depth+individualStepDepthDifference*i;
        if(DEPTHBUFFER[x][y] == 0){
            float pixel = (textureLine[round(i*oneStepThroughTexLine)].y * map.width) +textureLine[round(i*oneStepThroughTexLine)].x;
            uint32_t colour = map.pixels[pixel];
            window.setPixelColour(round(x), from.y, colour);
            DEPTHBUFFER[x][y] = 1/depth;         
        }else if(DEPTHBUFFER[x][y]>(1/depth)){
            float pixel = (textureLine[round(i*oneStepThroughTexLine)].y * map.width) +textureLine[round(i*oneStepThroughTexLine)].x;
            uint32_t colour = map.pixels[pixel];
            window.setPixelColour(round(x), from.y, colour);
            DEPTHBUFFER[x][y] = 1/depth;  
        }

    }
}

void drawTopHalfTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap map){
    vector<CanvasPoint> topSide1 = interpolateSideForTex(triangle.v0(), triangle.v1());
    vector<CanvasPoint> topSide2 = interpolateSideForTex(triangle.v0(), triangle.v2());
    vector<CanvasPoint> texTopside1 = interpolateSideForTex(triangle.v0().texturePoint, triangle.v1().texturePoint);
    vector<CanvasPoint> texTopside2 = interpolateSideForTex(triangle.v0().texturePoint, triangle.v2().texturePoint);
    for(float i=0; i<topSide1.size(); i++){
        for(float j=0; j<topSide2.size(); j++){
            if(topSide1[i].y == topSide2[j].y){
                float side1Proportion = i/topSide1.size();
                float side2Proportion = j/topSide2.size();
                float texturePosition1 = round(texTopside1.size()*side1Proportion);
                float texturePosition2 = round(texTopside2.size()*side2Proportion);
				drawTextureLine(texTopside1[texturePosition1],texTopside2[texturePosition2], topSide1[i], topSide2[j], map, window);
            }
        }
    }
}
void drawBottomHalfTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap map){
    vector<CanvasPoint> bottomSide1 = interpolateSideForTex(triangle.v0(), triangle.v2());
    vector<CanvasPoint> bottomSide2 = interpolateSideForTex(triangle.v1(), triangle.v2());  
    vector<CanvasPoint> texbottomSide1 = interpolateSideForTex(triangle.v0().texturePoint, triangle.v2().texturePoint);
    vector<CanvasPoint> texbottomSide2 = interpolateSideForTex(triangle.v1().texturePoint, triangle.v2().texturePoint); 
	for(float i=0; i<bottomSide1.size(); i++){
        for(float j=0; j<bottomSide2.size(); j++){
            if(bottomSide1[i].y == bottomSide2[j].y){
                float side1Proportion = i/bottomSide1.size();
                float side2Proportion = j/bottomSide2.size();
                float texturePosition1 = round(texbottomSide1.size()*side1Proportion);
                float texturePosition2 = round(texbottomSide2.size()*side2Proportion);
                drawTextureLine(texbottomSide1[texturePosition1],texbottomSide2[texturePosition2], bottomSide1[i], bottomSide2[j], map, window);
            }
        }
    }
}

void textured(DrawingWindow &window, CanvasTriangle triangle){
    TextureMap map = TextureMap("libs/sdw/texture.ppm");
    // triangle.v0().texturePoint = TexturePoint(195, 5);
    // triangle.v1().texturePoint = TexturePoint(395, 380);
    // triangle.v2().texturePoint = TexturePoint(65, 330);
	triangle.v0().texturePoint = TexturePoint(round(triangle.v0().texturePoint.x*map.width), round(map.height - triangle.v0().texturePoint.y*map.height));
    triangle.v1().texturePoint = TexturePoint(round(triangle.v1().texturePoint.x*map.width), round(map.height - triangle.v1().texturePoint.y*map.height));
    triangle.v2().texturePoint = TexturePoint(round(triangle.v2().texturePoint.x*map.width), round(map.height - triangle.v2().texturePoint.y*map.height));
	cout << triangle.v0().texturePoint << endl << triangle.v1().texturePoint << endl << triangle.v2().texturePoint << endl;
	triangle = sort(triangle);
	if(triangle.v0().y == triangle.v1().y){
		cout<< "drawing bottom half" << endl;
		drawBottomHalfTriangle(window, triangle, map);
	}
	else if(triangle.v1().y == triangle.v2().y){
		cout << "drawing top hafl" <<endl;
		drawTopHalfTriangle(window, triangle, map);
	}
	else{
		CanvasPoint newPoint = findNewPointTex(triangle);
		CanvasTriangle topTriangle = CanvasTriangle(triangle.v0(), triangle.v1(), newPoint);
		CanvasTriangle bottomTriangle = CanvasTriangle(newPoint, triangle.v1(), triangle.v2());
		drawTopHalfTriangle(window, topTriangle, map);
		drawBottomHalfTriangle(window, bottomTriangle, map);	
	}
}


void rasturisedDraw(DrawingWindow &window, vector<CanvasTriangle> canvasTriangles, vector<Colour> colours, vector<CanvasPoint> points, string choice){
	
	if(choice == "rasturising colour"){
		window.clearPixels();
		initialiseDepthBuffer();
		for(int i=0; i<canvasTriangles.size(); i++){
			Colour tempColour = colours[i];
			drawFilledTriangle(window, canvasTriangles[i], tempColour);
		}
		window.renderFrame();
	}
	else if(choice == "point"){
		window.clearPixels();
		initialiseDepthBuffer();
		for(int i=0; i<points.size(); i++){
			drawPoint(window, points[i]);
		}
		window.renderFrame();
	}
	else if(choice =="wireframe"){
		window.clearPixels();
		initialiseDepthBuffer();
		for(int i=0; i<canvasTriangles.size(); i++){
			Colour white = Colour(255,255,255);
			drawEmptyTriangle(window, canvasTriangles[i], white);
		}
		window.renderFrame();
	}
	else if(choice == "texturemap"){
		window.clearPixels();
		initialiseDepthBuffer();
		for(int i=0; i<canvasTriangles.size(); i++){
			Colour tempColour = colours[i];
			if(canvasTriangles[i].textured){
				canvasTriangles[i].v0().x = round(canvasTriangles[i].v0().x);
				canvasTriangles[i].v0().y = round(canvasTriangles[i].v0().y);
				canvasTriangles[i].v1().x = round(canvasTriangles[i].v1().x);
				canvasTriangles[i].v1().y = round(canvasTriangles[i].v1().y);
				canvasTriangles[i].v2().x = round(canvasTriangles[i].v2().x);
				canvasTriangles[i].v2().y = round(canvasTriangles[i].v2().y);
				textured(window, canvasTriangles[i]);
			}
			else{
				drawFilledTriangle(window, canvasTriangles[i], tempColour);
			}
        }
		window.renderFrame();
	}
	
}



//								************################# SECTION FOR RAY TRACE METHOD #################************

RayTriangleIntersection getClosestIntersection(glm::vec3 start, glm::vec3 direction, vector<ModelTriangle> triangles){
	vector<glm::vec4> viableIntersections;
	glm::vec4 closestIntersection;
	for(int i=0; i<triangles.size();i++){
		ModelTriangle triangle = triangles[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = start - triangle.vertices[0];
		glm::mat3 DEMatrix(-direction, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution[0];
		float u = possibleSolution[1];
		float v = possibleSolution[2];
		if(t>0 && u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0){
			viableIntersections.push_back(glm::vec4(t,u,v,i));
		}
	}
	for (int i=0; i<viableIntersections.size(); i++){
		if(i==0){
			closestIntersection = viableIntersections[0];
		}else{
			if(viableIntersections[i].x < closestIntersection.x){
				closestIntersection = viableIntersections[i];
			}
		}
	}
	RayTriangleIntersection intersectionTriangle;
	intersectionTriangle.intersectionPoint = start + (direction*closestIntersection[0]);
	intersectionTriangle.distanceFromCamera = closestIntersection[0];
	intersectionTriangle.intersectedTriangle = triangles[closestIntersection[3]];
	intersectionTriangle.triangleIndex = closestIntersection[3];
	intersectionTriangle.vertexNormals = triangles[closestIntersection[3]].vertexNormals;
	intersectionTriangle.u = closestIntersection[1];
	intersectionTriangle.v = closestIntersection[2];
	return intersectionTriangle;
}

bool doesShadowRayHitLight(RayTriangleIntersection originPoint, vector<ModelTriangle> modelTriangles, glm::vec3 lightposition){
	glm::vec3 direction = lightposition - originPoint.intersectionPoint;
	vector<ModelTriangle> newTriangles;
	for(int i = 0; i<modelTriangles.size();i++){
		if(i!=originPoint.triangleIndex){
			newTriangles.push_back(modelTriangles[i]);
		}
	}
	RayTriangleIntersection intersectionFromPointToLight = getClosestIntersection(originPoint.intersectionPoint, direction,newTriangles);
	if(intersectionFromPointToLight.distanceFromCamera != 0){
		if(abs(glm::l2Norm(intersectionFromPointToLight.intersectionPoint,originPoint.intersectionPoint))>abs(glm::l2Norm(LIGHTPOSITION,originPoint.intersectionPoint))){
			return true;
		}else{
			return false;
		}
	}else return true;
	
}

float howManyShadowRaysHitLight(RayTriangleIntersection originPoint, vector<ModelTriangle> modelTriangles){
	int numOfRaysHitLight = 0;
	for (int i=0; i<LIGHTPOSITIONS.size(); i++){
		bool hit = doesShadowRayHitLight(originPoint, modelTriangles, LIGHTPOSITIONS[i]);
		if(hit){
			numOfRaysHitLight++;
		}
	}
	return numOfRaysHitLight;
}

array<float, 3> getVertexBrightnesses(RayTriangleIntersection intersection, glm::vec3 cameraPosition){
	array<float, 3> returnArray;
	for(int i = 0; i<3; i++){
		glm::vec3 vectorFromLight = glm::normalize(LIGHTPOSITION - intersection.intersectionPoint);
		float aoiBrightness = glm::dot(intersection.vertexNormals[i], vectorFromLight);
		
		int n = 100;
		glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - LIGHTPOSITION);
		glm::vec3 normal = glm::normalize(intersection.vertexNormals[i]);
		glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*normal*(glm::dot(vectorOfLightToPoint, normal))));
		glm::vec3 vectorToCamera = glm::normalize((cameraPosition - intersection.intersectionPoint));
		float specBrightness = glm::dot(vectorToCamera, vectorOfReflection);
		if (specBrightness < 0){
			specBrightness = 0;
		}
		specBrightness = pow(specBrightness, n);
		
		float vertexBrightness = (specBrightness + aoiBrightness)/2;
		if (vertexBrightness<0){
			vertexBrightness = 0;
		}
		if (vertexBrightness>1){
			vertexBrightness = 1;
		}
		returnArray[i] = vertexBrightness;
	}
	return returnArray;
}

glm::vec3 getNewDirection(glm::vec3 rayDirection, glm::vec3 normal){
	glm::vec3 rayOfIncidence = glm::normalize(rayDirection);
	normal = glm::normalize(normal);
	glm::vec3 rayOfReflection = rayOfIncidence - (2*normal*(glm::dot(rayOfIncidence,normal)));
	rayOfReflection = glm::normalize(rayOfReflection);
	return rayOfReflection;
}

RayTriangleIntersection mirroredFunction(glm::vec3 rayDirection, RayTriangleIntersection intersection, vector<ModelTriangle> modelTriangles){
	glm::vec3 newDirection = getNewDirection(rayDirection, intersection.intersectedTriangle.normal);
	vector<ModelTriangle> newTriangles;
	for(int i = 0; i<modelTriangles.size();i++){
		if(i!=intersection.triangleIndex){
			newTriangles.push_back(modelTriangles[i]);
		}
	}
	RayTriangleIntersection newIntersection = getClosestIntersection(intersection.intersectionPoint, newDirection, newTriangles);
	if(intersection.triangleIndex<= newIntersection.triangleIndex){
		newIntersection.triangleIndex++;
	}
	return newIntersection;
}

void rayTraceDraw(DrawingWindow &window, glm::vec3 cameraPosition, vector<ModelTriangle> modelTriangles, vector<Colour> colours, float focalLength, string choice){
	if(choice == "rayTrace" || choice == "proximity" || choice == "angleOfIncidence" || choice == "specular" || choice == "ambient" || choice == "all" || choice == "gouraud" || choice == "phong" || choice == "soft") {
		window.clearPixels();
		float scaleDown = 400;
		uint32_t col;
		for(int y=0; y<HEIGHT; y++){
			for(int x=0; x<WIDTH; x++){
				float xraydirection = (x-WIDTH/2);
				float yraydirection = (y-HEIGHT/2);
				xraydirection = xraydirection / scaleDown;
				yraydirection = - yraydirection / scaleDown;
				glm::vec3 rayDirection = glm::vec3(xraydirection, yraydirection, cameraPosition.z - focalLength);
				rayDirection = rayDirection - cameraPosition;
				RayTriangleIntersection intersection = getClosestIntersection(cameraPosition, rayDirection, modelTriangles);	
				if(intersection.intersectedTriangle.mirrored){
					glm::vec3 rayOfIncidence = rayDirection;
					while(intersection.intersectedTriangle.mirrored){
						intersection = mirroredFunction(rayOfIncidence, intersection, modelTriangles);
						rayOfIncidence = getNewDirection(rayOfIncidence, intersection.intersectedTriangle.normal);
					}
				}
				if (intersection.distanceFromCamera != 0 ){
					if(choice == "rayTrace"){
						col = (255<<24) + (int(colours[intersection.triangleIndex].red) << 16) +(int(colours[intersection.triangleIndex].green) << 8) +int(colours[intersection.triangleIndex].blue);
					}
					else if (choice == "ambient"){
						if(doesShadowRayHitLight(intersection, modelTriangles, LIGHTPOSITION)){
							col = (255<<24) + (int(colours[intersection.triangleIndex].red) << 16) +(int(colours[intersection.triangleIndex].green) << 8) +int(colours[intersection.triangleIndex].blue);
						}
						else{
							float min = 0.2;
							col = (255<<24) + (int(colours[intersection.triangleIndex].red*min) << 16) +(int(colours[intersection.triangleIndex].green*min) << 8) +int(colours[intersection.triangleIndex].blue*min);
						}
					}
					else if (choice == "soft"){
						float numOfRaysHitLight = howManyShadowRaysHitLight(intersection, modelTriangles);
						double lightIntensity = numOfRaysHitLight/LIGHTPOSITIONS.size();
						if(lightIntensity<0.2){
							lightIntensity = 0.2;
						}
						col = (255<<24) + (int(colours[intersection.triangleIndex].red*lightIntensity) << 16) +(int(colours[intersection.triangleIndex].green*lightIntensity) << 8) +int(colours[intersection.triangleIndex].blue*lightIntensity);
					}
					else if (choice == "proximity"){
						glm::vec3 distanceToLightVec = LIGHTPOSITION - intersection.intersectionPoint;
						float distanceToLightVal = glm::length(distanceToLightVec);
						float mult = 2*1/(4*M_PI*pow(distanceToLightVal,2));
						if (mult > 1){
							mult = 1;
						}
						col = (255<<24) + (int((colours[intersection.triangleIndex].red)*mult) << 16) + (int((colours[intersection.triangleIndex].green)*mult) << 8) + int((colours[intersection.triangleIndex].blue)*mult);
					}
					else if (choice == "angleOfIncidence"){
						glm::vec3 directionToLight = LIGHTPOSITION - intersection.intersectionPoint;
						directionToLight = glm::normalize(directionToLight);
						glm::vec3 normal = glm::normalize(intersection.intersectedTriangle.normal);
						float dotProduct = glm::dot(directionToLight, normal);
						if (dotProduct < 0){
							dotProduct = 0;
						}
						col = (255<<24) + (int((colours[intersection.triangleIndex].red)*dotProduct) << 16) + (int((colours[intersection.triangleIndex].green)*dotProduct) << 8) + int((colours[intersection.triangleIndex].blue)*dotProduct);
					}
					else if (choice == "specular"){
						int n = 250;
						glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - LIGHTPOSITION);
						glm::vec3 normal = glm::normalize(intersection.intersectedTriangle.normal);
						glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*normal*(glm::dot(vectorOfLightToPoint, normal))));
						glm::vec3 vectorToCamera = glm::normalize((cameraPosition - intersection.intersectionPoint));
						float dotProduct = glm::dot(vectorToCamera, vectorOfReflection);
						if (dotProduct < 0){
							dotProduct = 0;
						}
						float multiplier = pow(dotProduct, n);
						col = (255<<24) + (int((colours[intersection.triangleIndex].red)*multiplier) << 16) + (int((colours[intersection.triangleIndex].green)*multiplier) << 8) + int((colours[intersection.triangleIndex].blue)*multiplier);
					}
					else if (choice == "all"){

						float numOfRaysHitLight = howManyShadowRaysHitLight(intersection, modelTriangles);
						double lightIntensity = numOfRaysHitLight/LIGHTPOSITIONS.size();

						glm::vec3 distanceToLightVec = LIGHTPOSITION - intersection.intersectionPoint;
						float distanceToLightVal = glm::length(distanceToLightVec);
						float proximityMult = 2*1/(4*M_PI*pow(distanceToLightVal,2));
						if (proximityMult>1){
							proximityMult = 1;
						}

						glm::vec3 directionToLight = LIGHTPOSITION - intersection.intersectionPoint;
						directionToLight = glm::normalize(directionToLight);
						glm::vec3 normal = glm::normalize(intersection.intersectedTriangle.normal);
						float angleOfIncidenceMult = glm::dot(directionToLight, normal);
						if (angleOfIncidenceMult < 0){
							angleOfIncidenceMult = 0;
						}

						int n = 32;
						glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - LIGHTPOSITION);
						glm::vec3 normal2 = glm::normalize(intersection.intersectedTriangle.normal);
						glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*normal2*(glm::dot(vectorOfLightToPoint, normal))));
						glm::vec3 vectorToCamera = glm::normalize((cameraPosition - intersection.intersectionPoint));
						float dotProduct = glm::dot(vectorToCamera, vectorOfReflection);
						if (dotProduct < 0){
							dotProduct = 0;
						}
						float specularMultiplier = pow(dotProduct, n);

						float pixelBrightness = ((proximityMult + angleOfIncidenceMult + specularMultiplier)/3) + 0.1;
						pixelBrightness = pixelBrightness*lightIntensity;
						if (pixelBrightness > 1){
							pixelBrightness = 1;
						}
						if(pixelBrightness<0.2){
							pixelBrightness = 0.2;
						}
						col = (255<<24) + (int((colours[intersection.triangleIndex].red)*pixelBrightness) << 16) + (int((colours[intersection.triangleIndex].green)*pixelBrightness) << 8) + int((colours[intersection.triangleIndex].blue)*pixelBrightness);
					}
					else if(choice == "gouraud"){
						array<float, 3> vertexBrightnesses = getVertexBrightnesses(intersection, cameraPosition);
						float u = intersection.u;
						float v = intersection.v;
						float w = 1-(u+v);
						float pixelBrightness = vertexBrightnesses[0]*w + vertexBrightnesses[1]*u + vertexBrightnesses[2]*v;
						pixelBrightness = pixelBrightness + 0.1;
						if (pixelBrightness>1){
							pixelBrightness = 1;
						}
						col = (255<<24) + (int(colours[intersection.triangleIndex].red*pixelBrightness) << 16) + (int(colours[intersection.triangleIndex].green*pixelBrightness) << 8) + (int(colours[intersection.triangleIndex].blue*pixelBrightness));
					}
					else if(choice == "phong"){
						array<glm::vec3, 3> vertexNormals = intersection.vertexNormals;
						float u = intersection.u;
						float v = intersection.v;
						float w = 1-(u+v); 
						glm::vec3 pixelNormal = glm::normalize(glm::normalize(vertexNormals[0])*w + glm::normalize(vertexNormals[1])*u + glm::normalize(vertexNormals[2])*v);

						glm::vec3 vectorFromLight = glm::normalize(LIGHTPOSITION - intersection.intersectionPoint);
						float aoiBrightness = glm::dot(pixelNormal, vectorFromLight);
						if(aoiBrightness<0){
							aoiBrightness = 0;
						}

						int n = 256;
						glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - LIGHTPOSITION);
						glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*pixelNormal*(glm::dot(vectorOfLightToPoint, pixelNormal))));
						glm::vec3 vectorToCamera = glm::normalize((cameraPosition - intersection.intersectionPoint));
						float specBrightness = glm::dot(vectorToCamera, vectorOfReflection);
						if (specBrightness < 0){
							specBrightness = 0;
						}
						specBrightness = pow(specBrightness, n);
						
						float pixelBrightness = (specBrightness + aoiBrightness);
						if (pixelBrightness<0){
							pixelBrightness = 0;
						}
						if (pixelBrightness>1){
							pixelBrightness = 1;
						}
						col = (255<<24) + (int(colours[intersection.triangleIndex].red*pixelBrightness) << 16) + (int(colours[intersection.triangleIndex].green*pixelBrightness) << 8) + (int(colours[intersection.triangleIndex].blue*pixelBrightness));
					}
				}
				else{
					col = 0;
				}	
				window.setPixelColour(x, y, col);
			}
		}
		window.renderFrame();
	}
	
}



//								************################# SECTION FOR MAIN FUNCTIONS #################************

struct InitReturn {
	vector<CanvasPoint> canvasPoints;
	vector<CanvasTriangle> canvasTriangles;
	vector<Colour> triangleColours;
	vector<ModelTriangle> modelTriangles;
	glm::vec3 cameraPosition;
	float focalLength;
} typedef InitReturn;


glm::mat3 makeMatrix(float angle, char selection){
	double pi = 2*acos(0.0);
	angle = angle *(pi/180);
	if (selection == 'y'){
		glm::vec3 col1 = glm::vec3(1.0, 0.0, 0.0);
		glm::vec3 col2 = glm::vec3(0.0, cosf(angle), sinf(angle));
		glm::vec3 col3 = glm::vec3(0.0, -sinf(angle), cosf(angle));
		glm::mat3 matrix = glm::mat3(col1,col2,col3);
		return matrix;
	}
	else if (selection == 'x'){
		glm::vec3 col1 = glm::vec3(cosf(angle), 0.0, -sinf(angle));
		glm::vec3 col2 = glm::vec3(0.0, 1.0, 0.0);
		glm::vec3 col3 = glm::vec3(sinf(angle), 0.0, cosf(angle));
		glm::mat3 matrix = glm::mat3(col1,col2,col3);
		return matrix;
	}
	else if (selection == 'w'){
		glm::vec3 col1 = glm::vec3(1.0, 0.0, 0.0);
		glm::vec3 col2 = glm::vec3(0.0, cosf(angle), sinf(angle));
		glm::vec3 col3 = glm::vec3(0.0, -sinf(angle), cosf(angle));
		glm::mat3 matrix = glm::mat3(col1,col2,col3);
		return matrix;
	}
	else if (selection == 'a'){
		glm::vec3 col1 = glm::vec3(cosf(angle), 0.0, -sinf(angle));
		glm::vec3 col2 = glm::vec3(0.0, 1.0, 0.0);
		glm::vec3 col3 = glm::vec3(sinf(angle), 0.0, cosf(angle));
		glm::mat3 matrix = glm::mat3(col1,col2,col3);
		return matrix;
	}
	else return glm::mat3(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	
}

string handleEvent(SDL_Event event, DrawingWindow &window, vector<CanvasPoint> points, vector<CanvasTriangle> canvasTriangles, vector<Colour> colours, string choice) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT){
			LIGHTPOSITION[2] = LIGHTPOSITION[2] + 0.1;
			initiliseLightPositions();
		}
		else if (event.key.keysym.sym == SDLK_RIGHT){
			LIGHTPOSITION[2] = LIGHTPOSITION[2] - 0.1;
			initiliseLightPositions();
		}
		else if (event.key.keysym.sym == SDLK_UP){
			LIGHTPOSITION[1] = LIGHTPOSITION[1] + 0.1;
			initiliseLightPositions();
		}
		else if (event.key.keysym.sym == SDLK_DOWN){
			LIGHTPOSITION[1] = LIGHTPOSITION[1] - 0.1;
			initiliseLightPositions();
		}
		// if (event.key.keysym.sym == SDLK_LEFT){
		// 	cout << "LEFT" << endl;
		// 	glm::vec3 change = glm::vec3(-0.10, 0.0, 0.0);
		// 	CAMERAPOSITION = CAMERAPOSITION + change;
		// }
		// else if (event.key.keysym.sym == SDLK_RIGHT){
		// 	cout << "RIGHT" << endl;
		// 	glm::vec3 change = glm::vec3(0.10, 0.0, 0.0);
		// 	CAMERAPOSITION = CAMERAPOSITION + change;
		// }
		// else if (event.key.keysym.sym == SDLK_UP){
		// 	cout << "UP" << endl;
		// 	glm::vec3 change = glm::vec3(0.0, 0.10, 0.0);
		// 	CAMERAPOSITION = CAMERAPOSITION + change;
		// }
		// else if (event.key.keysym.sym == SDLK_DOWN){
		// 	cout << "DOWN" << endl;
		// 	glm::vec3 change = glm::vec3(0.0, -0.10, 0.0);
		// 	CAMERAPOSITION = CAMERAPOSITION + change;
		// }
		else if (event.key.keysym.sym == SDLK_1){
			cout << "drawing points" << endl;
			choice = "point";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_2){
			cout << "drawing wireframe" << endl;
			choice ="wireframe";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_3){
			cout << "drawing using rasturising techinque" << endl;
			choice = "rasturising colour";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_4){
			cout << "drawing with texture map" << endl;
			choice = "texturemap";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_5){
			cout << "drawing using raytrace technique" << endl;
			choice = "rayTrace";
			return choice;		
		}
		else if (event.key.keysym.sym == SDLK_6){
			cout << "drawing with proximity" << endl;
			choice = "proximity";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_7){
			cout << "drawing with angle of incidence" << endl;
			choice = "angleOfIncidence";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_8){
			cout << "drawint using specular lighting" <<endl;
			choice = "specular";
			return choice;
		}
		else if (event.key.keysym.sym ==SDLK_9){
			cout <<"draw using ambient light" <<endl; 
			choice = "ambient";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_0){
			cout << "rendering with all lighting" << endl;
			choice = "all";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_MINUS){
			cout << "drawing using gouraud" << endl;
			choice  = "gouraud";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_EQUALS){
			cout << "drawing with phong shading" << endl;
			choice = "phong";
			return choice;
		}
		else if(event.key.keysym.sym == SDLK_z){
			cout << "drawing with soft shadows" << endl;
			choice = "soft";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_q){
			cout << "quitting..." << endl;
			Quit = true;
		}
	} 
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	choice = "nothing";
	cout << glm::to_string(CAMERAPOSITION) << endl;
	return choice;
}

glm::vec3 getNormalForATriangle(ModelTriangle triangle){
	glm::vec3 vertex1 = triangle.vertices[0];
	glm::vec3 vertex2 = triangle.vertices[1];
	glm::vec3 vertex3 = triangle.vertices[2];
	glm::vec3 side1 = vertex2 - vertex1;
	glm::vec3 side2 = vertex3 - vertex1;
	glm::vec3 crossProductVector = glm::cross(side1, side2);
	return crossProductVector;
}

vector<ModelTriangle> getNormalForAllTriangles(vector<ModelTriangle> modelTriangles){
	for(int i=0; i<modelTriangles.size(); i++){
		glm::vec3 normal = getNormalForATriangle(modelTriangles[i]);
		modelTriangles[i].normal = normal;
	}
	return modelTriangles;
}

//initialised all the canvas canvas points, triangles and colours ready for a render
InitReturn initialise(){
	vector<CanvasPoint> canvasPoints;
	vector<ModelTriangle> modelTriangles;
	vector<CanvasTriangle> canvasTriangles;
	vector<Colour> triangleColours;
	glm::vec3 cameraPosition = 	glm::vec3(0.0, 0.0, 2.5);
	float focalLength = 2.0;
	initialiseDepthBuffer();
	//gets all the triangles and points from the folder 
	modelTriangles = parser();
	modelTriangles = getNormalForAllTriangles(modelTriangles);
	for(int i=0; i<modelTriangles.size(); i++){
		for(int j=0; j<3; j++){
			CanvasPoint temp = getCanvasIntersectionPoint(cameraPosition, modelTriangles[i].vertices[j], focalLength);
            temp.texturePoint = modelTriangles[i].texturePoints[j];
			canvasPoints.push_back(temp);
		}
		triangleColours.push_back(modelTriangles[i].colour);
	}
	//puts the canvas points back into triangles but now with coords for the canvas
	for(float i=0; i<canvasPoints.size(); i=i+3){
		CanvasTriangle temp = CanvasTriangle(canvasPoints[i], canvasPoints[i+1], canvasPoints[i+2]);
		temp.textured = modelTriangles[i/3].textured;
		canvasTriangles.push_back(temp);
	}
	InitReturn final;
	final.canvasPoints = canvasPoints;
	final.canvasTriangles = canvasTriangles;
	final.triangleColours = triangleColours;
	final.modelTriangles = modelTriangles;
	final.cameraPosition = cameraPosition;
	final.focalLength = focalLength;
	return final;
}




void lookAt(){
	glm::vec3 forward = glm::normalize(CAMERAPOSITION);
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0),forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
	CAMERAORIENTATION =glm::mat3(right, up, forward);
}

void autoOrbit(DrawingWindow &window, string direction){
	for(int i=0; i<180;i++){
		cout << i << endl;
		InitReturn variables = initialise();
		window.clearPixels();
		initialiseDepthBuffer();
		glm::mat3 change;
		if(direction == "auto orbit down"){ change = makeMatrix(-2.0,'y');}
		if(direction == "auto orbit up"){ change = makeMatrix(2.0,'y');}
		if(direction == "auto orbit left"){ change = makeMatrix(2.0,'x');}
		if(direction == "auto orbit right"){ change = makeMatrix(-2.0,'x');}
		CAMERAPOSITION = CAMERAPOSITION*change;
		lookAt();
		for(int i=0; i<variables.canvasTriangles.size(); i++){
			Colour tempColour = variables.triangleColours[i];
			drawFilledTriangle(window, variables.canvasTriangles[i], tempColour);
		}
		window.renderFrame();
	}
}


int main(int argc, char *argv[]) {
	initiliseLightPositions();
	InitReturn variables = initialise();
	string choice;
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		if (window.pollForInputEvents(event)) {
			choice = handleEvent(event, window, variables.canvasPoints, variables.canvasTriangles, variables.triangleColours, choice);
			rasturisedDraw(window,variables.canvasTriangles,variables.triangleColours, variables.canvasPoints, choice);
			rayTraceDraw(window, variables.cameraPosition, variables.modelTriangles, variables.triangleColours, variables.focalLength, choice);
		} 
		if (Quit){break;}
	}
}

