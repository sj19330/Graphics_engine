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
glm::vec3 CAMERAPOSITION = 	glm::vec3(0.0, 0.0, 2.5);
glm::mat3 CAMERAORIENTATION = glm::mat3(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
glm::vec3 LIGHTPOSITION = glm::vec3(0.0, 0.4, 0.0);
vector<glm::vec3> LIGHTPOSITIONS;

float DEPTHBUFFER[WIDTH][HEIGHT];

//initalises depth buffer by setting all values to 0
void initialiseDepthBuffer(){
	for(int i=0; i<WIDTH; i++){
		for(int j=0; j<HEIGHT; j++){
			DEPTHBUFFER[i][j] = 0;
		}
	}
}
//initilises multiple light positions for soft shadows arround the centeral light position
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
	}
}

//								************################# SECTION FOR PARSER #################************
struct getNumbersFromLineReturn{
	glm::vec3 coords;
	glm::vec3 normals;
	glm::vec3 textures;
}typedef getNumbersFromLineReturn;

//reads the numbers from a line of a file
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
    }
	if (normals.size() == 3){
		returnValues.normals = glm::vec3(normals[0], normals[1], normals[2]);
	}else returnValues.normals = empty;
	return returnValues;
}

//reads the string from a line of a file (for colours)
string getWordFromLine(string str){
	stringstream line;
	string word;

	line << str;
	while(!line.eof()){
		line >> word;
	}
	return word;

}
//read in both mtl and obj file and return a vector of model triangles with their colours
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
	cfile.open("src/cornell-box.mtl");
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
	file.open("src/cornell-box.obj");
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

//								************################# SECTION FOR MAIN FUNCTIONS #################************

struct InitReturn {
	vector<CanvasPoint> canvasPoints;
	vector<CanvasTriangle> canvasTriangles;
	vector<Colour> triangleColours;
	vector<ModelTriangle> modelTriangles;
	float focalLength;
} typedef InitReturn;

//gets the canvas intersection points
CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float focalLength){
	float scaleUp = 400;
	glm::vec3 cameraToVertex = vertexPosition - CAMERAPOSITION;
	glm::vec3 adjusted = cameraToVertex*CAMERAORIENTATION;
	float depth = adjusted.z;
	float canvasX = (focalLength*((adjusted.x)/depth))*scaleUp +WIDTH/2;
	float canvasY = (focalLength*((adjusted.y)/depth))*scaleUp +HEIGHT/2;
	canvasX = WIDTH - canvasX;
	return CanvasPoint(canvasX, canvasY, depth);
}
// creates a matrix for changing orientation of the camera 
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

// handles the keypress events and can capture screenshots.
//to render pick the render you want and press the corresponding key
//when moving camera or orientation, press how you want to manipulate the camera then after press the render you want to see
//clicking just the camera manipulation buttons or move light buttons wont render anything untill a render button is pressed afterwards
//the sphere objects may not render correctly untill the light source his moved outside the sphere,
// press the left arrow key a few times then the render you want and it will show the sphere rendered properly
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
		else if(event.key.keysym.sym == SDLK_x){
			LIGHTPOSITION[0] = LIGHTPOSITION[0] - 0.1;
			initiliseLightPositions();
		}
		else if(event.key.keysym.sym == SDLK_c){
			LIGHTPOSITION[0] = LIGHTPOSITION[0] + 0.1;
			initiliseLightPositions();
		}
		else if (event.key.keysym.sym == SDLK_v){
			cout << "LEFT" << endl;
			glm::vec3 change = glm::vec3(-0.10, 0.0, 0.0);
			CAMERAPOSITION = CAMERAPOSITION + change;
		}
		else if (event.key.keysym.sym == SDLK_m){
			cout << "RIGHT" << endl;
			glm::vec3 change = glm::vec3(0.10, 0.0, 0.0);
			CAMERAPOSITION = CAMERAPOSITION + change;
		}
		else if (event.key.keysym.sym == SDLK_b){
			cout << "TOWARDS" << endl;
			glm::vec3 change = glm::vec3(0.0, 0.0, -0.1);
			CAMERAPOSITION = CAMERAPOSITION + change;
		}
		else if (event.key.keysym.sym == SDLK_n){
			cout << "BACK" << endl;
			glm::vec3 change = glm::vec3(0.0, 0.0, 0.1);
			CAMERAPOSITION = CAMERAPOSITION + change;
		}
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
		else if (event.key.keysym.sym == SDLK_v){
			cout << "forward" << endl;
			glm::vec3 change = glm::vec3(0.0, 0.0, -0.10);
			CAMERAPOSITION = CAMERAPOSITION + change;
		}
		else if (event.key.keysym.sym == SDLK_b){
			cout << "back" << endl;
			glm::vec3 change = glm::vec3(0.0, 0.0, 0.10);
			CAMERAPOSITION = CAMERAPOSITION + change;
		}
		else if (event.key.keysym.sym == SDLK_j){
			cout << "orbit left" << endl;
			glm::mat3 change = makeMatrix(5.0,'x');
			CAMERAPOSITION = CAMERAPOSITION*change;
		}
		else if (event.key.keysym.sym == SDLK_l){
			cout << "orbit right" << endl;
			glm::mat3 change = makeMatrix(-5.0,'x');
			CAMERAPOSITION = CAMERAPOSITION*change;
		}
		else if (event.key.keysym.sym == SDLK_i){
			cout << "orbit up" << endl;
			glm::mat3 change = makeMatrix(10.0,'y');
			CAMERAPOSITION = CAMERAPOSITION*change;
		}
		else if (event.key.keysym.sym == SDLK_k){
			cout << "orbit down" << endl;
			glm::mat3 change = makeMatrix(-5.0,'y');
			CAMERAPOSITION = CAMERAPOSITION*change;
		}
		else if(event.key.keysym.sym == SDLK_w){
			cout << "pan up" << endl;
			glm::mat3 rotationMatrix = makeMatrix(5.0,'w');
			CAMERAORIENTATION = rotationMatrix*CAMERAORIENTATION;
		}
		else if(event.key.keysym.sym == SDLK_s){
			cout << "pan down" << endl;
			glm::mat3 rotationMatrix = makeMatrix(-10.0,'w');
			CAMERAORIENTATION = rotationMatrix*CAMERAORIENTATION;
		}
		else if (event.key.keysym.sym == SDLK_a){
			cout <<"pan left" << endl;
			glm::mat3 rotationMatrix = makeMatrix(5.0,'a');
			CAMERAORIENTATION = rotationMatrix*CAMERAORIENTATION;
		}
		else if (event.key.keysym.sym == SDLK_d){
			cout <<"pan right" << endl;
			glm::mat3 rotationMatrix = makeMatrix(-5.0,'a');
			CAMERAORIENTATION = rotationMatrix*CAMERAORIENTATION;
		}
		else if (event.key.keysym.sym == SDLK_g){
			cout <<"auto orbit down" << endl;
			choice = "auto orbit down";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_t){
			cout <<"auto orbit up" << endl;
			choice = "auto orbit up";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_h){
			cout <<"auto orbit right" << endl;
			choice = "auto orbit right";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_f){
			cout <<"auto orbit left" << endl;
			choice = "auto orbit left";
			return choice;
		}
		else if (event.key.keysym.sym == SDLK_q){
			cout << "quitting..." << endl;
			Quit = true;
		}
	} 
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		string fileName = "ouput.ppm";
		window.savePPM(fileName);
	}
	choice = "nothing";
	return choice;
}

//returns the vector of the normal of a model triangle
glm::vec3 getNormalForATriangle(ModelTriangle triangle){
	glm::vec3 vertex1 = triangle.vertices[0];
	glm::vec3 vertex2 = triangle.vertices[1];
	glm::vec3 vertex3 = triangle.vertices[2];
	glm::vec3 side1 = vertex2 - vertex1;
	glm::vec3 side2 = vertex3 - vertex1;
	glm::vec3 crossProductVector = glm::cross(side1, side2);
	return crossProductVector;
}

//loops through all model triangles, finds their normals and addes it to the model
//triangle normal section before returning all the model triangles back
vector<ModelTriangle> getNormalForAllTriangles(vector<ModelTriangle> modelTriangles){
	for(int i=0; i<modelTriangles.size(); i++){
		glm::vec3 normal = getNormalForATriangle(modelTriangles[i]);
		modelTriangles[i].normal = normal;
	}
	return modelTriangles;
}

//initialises all the canvas canvas points, triangles and colours ready for a render
InitReturn initialise(){
	vector<CanvasPoint> canvasPoints;
	vector<ModelTriangle> modelTriangles;
	vector<CanvasTriangle> canvasTriangles;
	vector<Colour> triangleColours;
	float focalLength = 2.0;
	initialiseDepthBuffer();
	//gets all the triangles and points from the folder 
	modelTriangles = parser();
	modelTriangles = getNormalForAllTriangles(modelTriangles);
	for(int i=0; i<modelTriangles.size(); i++){
		for(int j=0; j<3; j++){
			CanvasPoint temp = getCanvasIntersectionPoint(modelTriangles[i].vertices[j], focalLength);
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
	final.focalLength = focalLength;
	return final;
}








//								************################# SECTION FOR RASTURISING METHOD #################************



//draws a white point on screen for a canvaspoint
void drawPoint(DrawingWindow &window, CanvasPoint point) {
	uint32_t white = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	window.setPixelColour(point.x, point.y, white);
}

//draws a coloured line from a canvas point to another
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
// draws the outline of a canvas triangle
void drawEmptyTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour col){
	drawLineConsideringDepth(window, triangle.v1(), triangle.v0(), col);
	drawLineConsideringDepth(window, triangle.v1(), triangle.v2(), col);
	drawLineConsideringDepth(window, triangle.v2(), triangle.v0(), col);
}
//sorts a canvas triangle so the point with the middle y valuse is put as v1
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
//finds the newpoint for a canvas triangle on the v0,v2 side, in line with the y value of v1
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
//interpolates a line from one point to another, returning all the canvas points inbetween 
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
//interpolates the two sides of canvas triangle, then fills the triangle in with colour
void interpolateSidesAndFill(CanvasTriangle tri, DrawingWindow &window, Colour col, int side){
	vector<CanvasPoint> side1 = interpolateLine(tri.v0(), tri.v1(), side);
	vector<CanvasPoint> side2 = interpolateLine(tri.v0(), tri.v2(), side);
	for(int i = 0; i< side1.size(); i++){
		drawLineConsideringDepth(window, side1[i], side2[i], col);
	}
	drawEmptyTriangle(window, tri, col);
}
//draws a canvas triangle, filling it with a specified colour
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


//interpolates from a start canvaspoint to an end canvaspoint for textured triangles
vector<CanvasPoint> interpolateSideForTex(CanvasPoint start, CanvasPoint end){
    vector<CanvasPoint> v;
    float numberOfValuesy = abs(round(end.y) - round(start.y));
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
// same as above function but for texture points
vector<CanvasPoint> interpolateSideForTex(TexturePoint start, TexturePoint end){
    vector<CanvasPoint> v;
    float numberOfTexValuesy = abs(round(end.y) - round(start.y));
    float numberOfTexValuesx = abs(round(end.x) - round(start.x));
    float numberOfTexValues = fmax(numberOfTexValuesy, numberOfTexValuesx);
    float texStepSizex = (end.x - start.x)/numberOfTexValues;
    float texStepSizey = (end.y - start.y)/numberOfTexValues;
    for (int i=0; i<numberOfTexValues; i++){
        CanvasPoint temp = CanvasPoint(round(start.x + i*texStepSizex), round(start.y + i *texStepSizey));
        v.push_back(temp);
    }
    return v;
}
//finds the new point for a textured triangle
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

//draws a line using pixels from a texture map from a start to an end point
void drawTextureLine(CanvasPoint texStart, CanvasPoint texEnd, CanvasPoint from, CanvasPoint to, TextureMap map, DrawingWindow &window, string topOrBottom){
	float numberOfSteps = abs(to.x - from.x);
    float depthDiff = to.depth - from.depth;
	float individualStepDepthDifference = depthDiff/numberOfSteps;
    vector<CanvasPoint> textureLine = interpolateSideForTex(texStart, texEnd);
    for(float i=0.0; i<numberOfSteps; i++){
        float oneStepThroughTexLine = (textureLine.size()-1)/numberOfSteps;
        int x = 0;
		if (topOrBottom == "bottom"){
			x = from.x - i;
		}
		if(topOrBottom == "top"){
			x = from.x + i;
		}
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
//draws the top half of a textured triangle
void drawTopHalfTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap map){
    vector<CanvasPoint> topSide1 = interpolateSideForTex(triangle.v0(), triangle.v2());
    vector<CanvasPoint> topSide2 = interpolateSideForTex(triangle.v0(), triangle.v1());
    vector<CanvasPoint> texTopside1 = interpolateSideForTex(triangle.v0().texturePoint, triangle.v2().texturePoint);
    vector<CanvasPoint> texTopside2 = interpolateSideForTex(triangle.v0().texturePoint, triangle.v1().texturePoint);
    for(float i=0; i<topSide1.size(); i++){
        for(float j=0; j<topSide2.size(); j++){
            if(topSide1[i].y == topSide2[j].y){
                float side1Proportion = i/topSide1.size();
                float side2Proportion = j/topSide2.size();
                float texturePosition1 = round(texTopside1.size()*side1Proportion);
                float texturePosition2 = round(texTopside2.size()*side2Proportion);
				drawTextureLine(texTopside1[texturePosition1],texTopside2[texturePosition2], topSide1[i], topSide2[j], map, window, "top");
            }
        }
    }
}
//draws the bottom half of a textured triangle
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
                drawTextureLine(texbottomSide1[texturePosition1],texbottomSide2[texturePosition2], bottomSide1[i], bottomSide2[j], map, window, "bottom");
            }
        }
    }
}
// draws a textured triangle
void textured(DrawingWindow &window, CanvasTriangle triangle){
    TextureMap map = TextureMap("libs/sdw/texture.ppm");
	triangle.v0().texturePoint = TexturePoint(round(triangle.v0().texturePoint.x*map.width), round(map.height - triangle.v0().texturePoint.y*map.height));
    triangle.v1().texturePoint = TexturePoint(round(triangle.v1().texturePoint.x*map.width), round(map.height - triangle.v1().texturePoint.y*map.height));
    triangle.v2().texturePoint = TexturePoint(round(triangle.v2().texturePoint.x*map.width), round(map.height - triangle.v2().texturePoint.y*map.height));
	triangle = sort(triangle);
	CanvasPoint newPoint = findNewPointTex(triangle);
	CanvasTriangle topTriangle = CanvasTriangle(triangle.v0(), triangle.v1(), newPoint);
	CanvasTriangle bottomTriangle = CanvasTriangle(newPoint, triangle.v1(), triangle.v2());
	drawTopHalfTriangle(window, topTriangle, map);
	drawBottomHalfTriangle(window, bottomTriangle, map);	
}
// when called, lookat points the camera towards the the centre of the model
void lookAt(){
	glm::vec3 forward = glm::normalize(CAMERAPOSITION);
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0),forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
	CAMERAORIENTATION =glm::mat3(right, up, forward);
}
// drawing function for all methods that use the rasturised techinque
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
	else if (choice == "auto orbit down"){
		for (int i=0; i<180; i++){
			InitReturn variables = initialise();
			window.clearPixels();
			glm::mat3 change = makeMatrix(-2.0, 'y');
			CAMERAPOSITION = CAMERAPOSITION * change;
			lookAt();
			for(int i=0; i<variables.canvasTriangles.size(); i++){
				Colour tempColour = variables.triangleColours[i];
				drawEmptyTriangle(window, variables.canvasTriangles[i], tempColour);
			}
			window.renderFrame();
		}
	}
	else if (choice == "auto orbit up"){
		for (int i=0; i<180; i++){
			InitReturn variables = initialise();
			window.clearPixels();
			glm::mat3 change = makeMatrix(2.0, 'y');
			CAMERAPOSITION = CAMERAPOSITION * change;
			lookAt();
			for(int i=0; i<variables.canvasTriangles.size(); i++){
				Colour tempColour = variables.triangleColours[i];
				drawFilledTriangle(window, variables.canvasTriangles[i], tempColour);
			}
			window.renderFrame();
		}
	}
	else if (choice == "auto orbit left"){
		for (int i=0; i<180; i++){
			InitReturn variables = initialise();
			window.clearPixels();
			glm::mat3 change = makeMatrix(2.0, 'x');
			CAMERAPOSITION = CAMERAPOSITION * change;
			lookAt();
			for(int i=0; i<variables.canvasTriangles.size(); i++){
				Colour tempColour = variables.triangleColours[i];
				drawFilledTriangle(window, variables.canvasTriangles[i], tempColour);
			}
			window.renderFrame();
		}
	}
	else if (choice == "auto orbit right"){
		for (int i=0; i<180; i++){
			InitReturn variables = initialise();
			window.clearPixels();
			glm::mat3 change = makeMatrix(-2.0, 'x');
			CAMERAPOSITION = CAMERAPOSITION * change;
			lookAt();
			for(int i=0; i<variables.canvasTriangles.size(); i++){
				Colour tempColour = variables.triangleColours[i];
				drawFilledTriangle(window, variables.canvasTriangles[i], tempColour);
			}
			window.renderFrame();
		}
	}
}



//								************################# SECTION FOR RAY TRACE METHOD #################************
//gets the closest intersection of a model triangle from a given start point and a given direction
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
// returns wether a point can see the main light point 
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
//returns how many light points out of the multiple light points a certain point can see
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

//finds the brighness of each vertex of a triangle of an intersection using angle of incidence and specular lighting techniques
//this is used for gauroud shading
array<float, 3> getVertexBrightnesses(RayTriangleIntersection intersection){
	array<float, 3> returnArray;
	for(int i = 0; i<3; i++){
		glm::vec3 vectorFromLight = glm::normalize(LIGHTPOSITION - intersection.intersectionPoint);
		float aoiBrightness = glm::dot(intersection.vertexNormals[i], vectorFromLight);
		
		int n = 100;
		glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - LIGHTPOSITION);
		glm::vec3 normal = glm::normalize(intersection.vertexNormals[i]);
		glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*normal*(glm::dot(vectorOfLightToPoint, normal))));
		glm::vec3 vectorToCamera = glm::normalize((CAMERAPOSITION - intersection.intersectionPoint));
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

//finds the new direction of a ray that bounces off a mirrored surface
glm::vec3 getNewDirection(glm::vec3 rayDirection, glm::vec3 normal){
	glm::vec3 rayOfIncidence = glm::normalize(rayDirection);
	normal = glm::normalize(normal);
	glm::vec3 rayOfReflection = rayOfIncidence - (2*normal*(glm::dot(rayOfIncidence,normal)));
	rayOfReflection = glm::normalize(rayOfReflection);
	return rayOfReflection;
}

//funciton for mirrored surfaces. it finds the direction of the ray of reflection then returns 
//the closest intersection from the mirrored surface in the new direction
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

//draw funciton for all ray trace methods
void rayTraceDraw(DrawingWindow &window, vector<ModelTriangle> modelTriangles, vector<Colour> colours, float focalLength, string choice){
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
				glm::vec3 rayDirection = glm::vec3(xraydirection, yraydirection, CAMERAPOSITION.z - focalLength);
				rayDirection = rayDirection - CAMERAPOSITION;
				RayTriangleIntersection intersection = getClosestIntersection(CAMERAPOSITION, rayDirection, modelTriangles);	
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
						glm::vec3 vectorToCamera = glm::normalize((CAMERAPOSITION - intersection.intersectionPoint));
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

						int n = 64;
						glm::vec3 vectorOfLightToPoint = glm::normalize(intersection.intersectionPoint - LIGHTPOSITION);
						glm::vec3 normal2 = glm::normalize(intersection.intersectedTriangle.normal);
						glm::vec3 vectorOfReflection = glm::normalize((vectorOfLightToPoint - 2*normal2*(glm::dot(vectorOfLightToPoint, normal))));
						glm::vec3 vectorToCamera = glm::normalize((CAMERAPOSITION - intersection.intersectionPoint));
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
						array<float, 3> vertexBrightnesses = getVertexBrightnesses(intersection);
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
						glm::vec3 vectorToCamera = glm::normalize((CAMERAPOSITION - intersection.intersectionPoint));
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



//main function
int main(int argc, char *argv[]) {
	initiliseLightPositions();
	string choice;
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		InitReturn variables = initialise();
		if (window.pollForInputEvents(event)) {
			choice = handleEvent(event, window, variables.canvasPoints, variables.canvasTriangles, variables.triangleColours, choice);
			rasturisedDraw(window,variables.canvasTriangles,variables.triangleColours, variables.canvasPoints, choice);
			rayTraceDraw(window, variables.modelTriangles, variables.triangleColours, variables.focalLength, choice);
		} 
		if (Quit){break;}
	}
}