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

#define WIDTH 320*3
#define HEIGHT 240*3
using namespace std;
bool Quit = false;
glm::vec3 LIGHTPOSITION = glm::vec3(0.0, 0.0, 0.0);
int test = 1;


//								************################# SECTION FOR PARSER #################************
glm::vec3 getNumbersFromLine(string str){
	stringstream line;
	string word;
	float coord;
	vector<float> coords;
	
	line << str;
	while(!line.eof()){
		line >> word;
		if(stringstream(word) >> coord){
			coords.push_back(coord);
		}
	}
	return glm::vec3(coords[0],coords[1],coords[2]);
}
string getColourFromLine(string str){
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
	vector<glm::vec3> empty;
	vector<glm::vec3> flines;
	vector<ModelTriangle> triangles;
	ifstream file;
	ifstream cfile;
	unordered_map<string, glm::vec3> colours;
	string colourForThisObject = " ";

	//open read and creat a hashmap of mtl file contents then close
	cfile.open("src/cornell-box.mtl");
	if(!cfile.is_open()){
		cout << "error opening mtl file" << endl;
	}
	else{
		while(cfile.good()){
			getline(cfile, cline);
			if(cline[0]=='n'){
				string tempColour = getColourFromLine(cline);
				getline(cfile, cline);
				glm::vec3 rgb = getNumbersFromLine(cline);
				colours[tempColour] = rgb;
			}
		}
	}
	cfile.close();
	//open and read the contents of the obj file then close
	file.open("src/cornell-box.obj");
	if (!file.is_open()){
		cout << "error opening obj file" << endl;
	}
	else{
		while(file.good()){
			getline(file, line);
			if(line[0] == 'u'){
				colourForThisObject = getColourFromLine(line);
			}
			else if(line[0] == 'v' && line[1] == ' '){
				glm::vec3 vertex = getNumbersFromLine(line);
				vertexes.push_back(vertex);
			}
			else if(line[0] == 'v' && line[1] == 'n'){
				glm::vec3 vertexNormal = getNumbersFromLine(line);
				vertexNormals.push_back(vertexNormal);
			}
			else if(line[0] == 'f'){
				glm::vec3 fline = getNumbersFromLine(line);
				Colour col = Colour(colours[colourForThisObject][0]*255, colours[colourForThisObject][1]*255, colours[colourForThisObject][2]*255);
				ModelTriangle temp = ModelTriangle(vertexes[((fline[0])-1)]*scale,vertexes[((fline[1])-1)]*scale,vertexes[((fline[2])-1)]*scale, col);
				if (vertexNormals != empty){
					temp.vertexNormals = {vertexNormals[(fline[0])-1], vertexNormals[(fline[1])-1], vertexNormals[(fline[2])-1]};
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
	float scaleUp = 400;
	//check if i should be returning depth of vertex position z PLUS camera position
	float depth = vertexPosition.z - cameraPosition.z;
	float canvasX = (focalLength*((vertexPosition.x)/depth))*scaleUp +WIDTH/2;
	float canvasY = (focalLength*((vertexPosition.y)/depth))*scaleUp +HEIGHT/2;
	canvasX = WIDTH - canvasX;
	depth = vertexPosition.z - cameraPosition.z;
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




bool doesShadowRayHitLight(RayTriangleIntersection originPoint, vector<ModelTriangle> modelTriangles){
	glm::vec3 direction = LIGHTPOSITION - originPoint.intersectionPoint;
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


void rayTraceDraw(DrawingWindow &window, glm::vec3 cameraPosition, vector<ModelTriangle> modelTriangles, vector<Colour> colours, float focalLength, string choice){
	if(choice == "rayTrace" || choice == "proximity" || choice == "hardShadow" || choice == "angleOfIncidence" || choice == "specular" || choice == "ambient" || choice == "all" || choice == "gouraud") {
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
				if (intersection.distanceFromCamera != 0 ){
					if(choice == "rayTrace"){
						col = (255<<24) + (int(colours[intersection.triangleIndex].red) << 16) +(int(colours[intersection.triangleIndex].green) << 8) +int(colours[intersection.triangleIndex].blue);
					}
					else if(choice == "hardShadow"){
						if(doesShadowRayHitLight(intersection, modelTriangles)){
							col = (255<<24) + (int(colours[intersection.triangleIndex].red) << 16) +(int(colours[intersection.triangleIndex].green) << 8) +int(colours[intersection.triangleIndex].blue);
						}
						else {
							col = 0;
						}
					}
					else if (choice == "ambient"){
						if(doesShadowRayHitLight(intersection, modelTriangles)){
							col = (255<<24) + (int(colours[intersection.triangleIndex].red) << 16) +(int(colours[intersection.triangleIndex].green) << 8) +int(colours[intersection.triangleIndex].blue);
						}
						else {
							float min = 0.2;
							col = (255<<24) + (int(colours[intersection.triangleIndex].red*min) << 16) +(int(colours[intersection.triangleIndex].green*min) << 8) +int(colours[intersection.triangleIndex].blue*min);
						}
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
						if(doesShadowRayHitLight(intersection, modelTriangles)){

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

							float pixelBrightness = (proximityMult + angleOfIncidenceMult + specularMultiplier)/3;
							col = (255<<24) + (int((colours[intersection.triangleIndex].red)*pixelBrightness) << 16) + (int((colours[intersection.triangleIndex].green)*pixelBrightness) << 8) + int((colours[intersection.triangleIndex].blue)*pixelBrightness);
						
						}
						else {
							float min = 0.2;
							col = (255<<24) + (int(colours[intersection.triangleIndex].red*min) << 16) +(int(colours[intersection.triangleIndex].green*min) << 8) +int(colours[intersection.triangleIndex].blue*min);
						}

					}
					else if(choice == "gouraud"){
						array<float, 3> vertexBrightnesses = getVertexBrightnesses(intersection, cameraPosition);
						float u = intersection.u;
						float v = intersection.v;
						float w = 1-(u+v);
						float pixelBrightness = vertexBrightnesses[0]*w + vertexBrightnesses[1]*u + vertexBrightnesses[2]*v;
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

string handleEvent(SDL_Event event, DrawingWindow &window, vector<CanvasPoint> points, vector<CanvasTriangle> canvasTriangles, vector<Colour> colours, string choice) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT){
			LIGHTPOSITION[2] = LIGHTPOSITION[2] + 0.1;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT){
			LIGHTPOSITION[2] = LIGHTPOSITION[2] - 0.1;
		}
		else if (event.key.keysym.sym == SDLK_UP){
			LIGHTPOSITION[1] = LIGHTPOSITION[1] + 0.1;
		}
		else if (event.key.keysym.sym == SDLK_DOWN){
			LIGHTPOSITION[1] = LIGHTPOSITION[1] - 0.1;
		}
		else if (event.key.keysym.sym == SDLK_p){
			cout << "drawing points" << endl;
			choice = "point";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_w){
			cout << "drawing wireframe" << endl;
			choice ="wireframe";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_t){
			cout << "drawing using raytrace technique" << endl;
			choice = "rayTrace";
			return choice;			
		}
		else if (event.key.keysym.sym == SDLK_h){
			cout << "drawing with hard shadows" << endl;
			choice = "hardShadow";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_x){
			cout << "drawing with proximity" << endl;
			choice = "proximity";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_i){
			cout << "drawing with angle of incidence" << endl;
			choice = "angleOfIncidence";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_s){
			cout << "drawint using specular lighting" <<endl;
			choice = "specular";
			return choice;
		}
		else if (event.key.keysym.sym ==SDLK_a){
			cout <<"draw using ambient light" <<endl; 
			choice = "ambient";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_r){
			cout << "drawing using rasturising techinque" << endl;
			choice = "rasturising colour";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_z){
			cout << "rendering with all lighting" << endl;
			choice = "all";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_g){
			cout << "drawing using gouraud" << endl;
			choice  = "gouraud";
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
	glm::vec3 cameraPosition = 	glm::vec3(0.0, 0.0, 1.0);
	float focalLength = 1.0;
	initialiseDepthBuffer();
	//gets all the triangles and points from the folder 
	modelTriangles = parser();
	modelTriangles = getNormalForAllTriangles(modelTriangles);
	for(int i=0; i<modelTriangles.size(); i++){
		for(int j=0; j<3; j++){
			CanvasPoint temp = getCanvasIntersectionPoint(cameraPosition, modelTriangles[i].vertices[j], focalLength);
			canvasPoints.push_back(temp);
		}
		triangleColours.push_back(modelTriangles[i].colour);
	}
	//puts the canvas points back into triangles but now with coords for the canvas
	for(int i=0; i<canvasPoints.size(); i=i+3){
		CanvasTriangle temp = CanvasTriangle(canvasPoints[i], canvasPoints[i+1], canvasPoints[i+2]);
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





int main(int argc, char *argv[]) {
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


