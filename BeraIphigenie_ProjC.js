//3456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_
// (JT: why the numbers? counts columns, helps me keep 80-char-wide listings)
//
// TABS set to 2.
//
// ORIGINAL SOURCE:
// RotatingTranslatedTriangle.js (c) 2012 matsuda
// HIGHLY MODIFIED to make:
//
// JT_MultiShader.js  for EECS 351-1, 
//									Northwestern Univ. Jack Tumblin

// Jack Tumblin's Project C -- step by step.

/* Show how to use 3 separate VBOs with different verts, attributes & uniforms. 
-------------------------------------------------------------------------------
	Create a 'VBObox' object/class/prototype & library to collect, hold & use all 
	data and functions we need to render a set of vertices kept in one Vertex 
	Buffer Object (VBO) on-screen, including:
	--All source code for all Vertex Shader(s) and Fragment shader(s) we may use 
		to render the vertices stored in this VBO;
	--all variables needed to select and access this object's VBO, shaders, 
		uniforms, attributes, samplers, texture buffers, and any misc. items. 
	--all variables that hold values (uniforms, vertex arrays, element arrays) we 
	  will transfer to the GPU to enable it to render the vertices in our VBO.
	--all user functions: init(), draw(), adjust(), reload(), empty(), restore().
	Put all of it into 'JT_VBObox-Lib.js', a separate library file.

USAGE:
------
1) If your program needs another shader program, make another VBObox object:
 (e.g. an easy vertex & fragment shader program for drawing a ground-plane grid; 
 a fancier shader program for drawing Gouraud-shaded, Phong-lit surfaces, 
 another shader program for drawing Phong-shaded, Phong-lit surfaces, and
 a shader program for multi-textured bump-mapped Phong-shaded & lit surfaces...)
 
 HOW:
 a) COPY CODE: create a new VBObox object by renaming a copy of an existing 
 VBObox object already given to you in the VBObox-Lib.js file. 
 (e.g. copy VBObox1 code to make a VBObox3 object).

 b) CREATE YOUR NEW, GLOBAL VBObox object.  
 For simplicity, make it a global variable. As you only have ONE of these 
 objects, its global scope is unlikely to cause confusions/errors, and you can
 avoid its too-frequent use as a function argument.
 (e.g. above main(), write:    var phongBox = new VBObox3();  )

 c) INITIALIZE: in your JS progam's main() function, initialize your new VBObox;
 (e.g. inside main(), write:  phongBox.init(); )

 d) DRAW: in the JS function that performs all your webGL-drawing tasks, draw
 your new VBObox's contents on-screen. 
 (NOTE: as it's a COPY of an earlier VBObox, your new VBObox's on-screen results
  should duplicate the initial drawing made by the VBObox you copied.  
  If that earlier drawing begins with the exact same initial position and makes 
  the exact same animated moves, then it will hide your new VBObox's drawings!
  --THUS-- be sure to comment out the earlier VBObox's draw() function call  
  to see the draw() result of your new VBObox on-screen).
  (e.g. inside drawAll(), add this:  
      phongBox.switchToMe();
      phongBox.draw();            )

 e) ADJUST: Inside the JS function that animates your webGL drawing by adjusting
 uniforms (updates to ModelMatrix, etc) call the 'adjust' function for each of your
VBOboxes.  Move all the uniform-adjusting operations from that JS function into the
'adjust()' functions for each VBObox. 

2) Customize the VBObox contents; add vertices, add attributes, add uniforms.
 ==============================================================================*/


// Global Variables  
//   (These are almost always a BAD IDEA, but here they eliminate lots of
//    tedious function arguments. 
//    Later, collect them into just a few global, well-organized objects!)
// ============================================================================
// for WebGL usage:--------------------
var gl;													// WebGL rendering context -- the 'webGL' object
																// in JavaScript with all its member fcns & data
var g_canvasID;									// HTML-5 'canvas' element ID#

// For multiple VBOs & Shaders:-----------------
worldBox = new VBObox0();		  // Holds VBO & shaders for 3D 'world' ground-plane grid, etc;
gouraudBox = new VBObox1();		  // "  "  for first set of custom-shaded 3D parts
phongBox = new VBObox2();     // "  "  for second set of custom-shaded 3D parts

// For animation:---------------------
var g_lastMS = Date.now();			// Timestamp (in milliseconds) for our 
                                // most-recently-drawn WebGL screen contents.  
                                // Set & used by moveAll() fcn to update all
                                // time-varying params for our webGL drawings.
  // All time-dependent params (you can add more!)

var materialNames = [
    "MATL_RED_PLASTIC",
    "MATL_GRN_PLASTIC",
    "MATL_BLU_PLASTIC",
    "MATL_BLACK_PLASTIC",
    "MATL_BLACK_RUBBER",
    "MATL_BRASS",
    "MATL_BRONZE_DULL",
    "MATL_BRONZE_SHINY",
    "MATL_CHROME",
    "MATL_COPPER_DULL",
    "MATL_COPPER_SHINY",
    "MATL_GOLD_DULL",
    "MATL_GOLD_SHINY",
    "MATL_PEWTER",
    "MATL_SILVER_DULL",
    "MATL_SILVER_SHINY",
    "MATL_EMERALD",
    "MATL_JADE",
    "MATL_OBSIDIAN",
    "MATL_PEARL",
    "MATL_RUBY",
    "MATL_TURQUOISE",
    "MATL_DEFAULT"
];
var selectedMaterial1;
var selectedInd1 = 0;

var selectedMaterial2;
var selectedInd2 = 0;

var g_currMatl = 1;
var matl1 = new Material();
var matl2 = new Material();


// Get the dropdown element
var materialDropdown1 = document.getElementById("materialDropdown1");
var materialDropdown2 = document.getElementById("materialDropdown2");

  // Populate dropdown with material names
  materialNames.forEach(function(materialName) {
    var option = document.createElement("option");
    option.text = materialName;
    materialDropdown1.add(option);
  });

   // Populate dropdown with material names
   materialNames.forEach(function(materialName) {
    var option = document.createElement("option");
    option.text = materialName;
    materialDropdown2.add(option);
  });


  // Add event listener for change
  materialDropdown1.addEventListener("change", function() {
    selectedMaterial1 = materialDropdown1.options[materialDropdown1.selectedIndex].text;
    selectedInd1 = materialDropdown1.selectedIndex;
    //document.getElementById("selectedMaterial").innerText = "Selected material: " + selectedMaterial+selectedInd.toString();
  });


  materialDropdown2.addEventListener("change", function() {
    selectedMaterial2 = materialDropdown2.options[materialDropdown2.selectedIndex].text;
    selectedInd2 = materialDropdown2.selectedIndex;
    //document.getElementById("selectedMaterial").innerText = "Selected material: " + selectedMaterial+selectedInd.toString();
  });

var g_angleNow0  =  0.0; 			  // Current rotation angle, in degrees.
var g_angleRate0 = 45.0;				// Rotation angle rate, in degrees/second.
                                //---------------
var g_angleNow1  = 100.0;       // current angle, in degrees
var g_angleRate1 =  95.0;        // rotation angle rate, degrees/sec
var g_angleMax1  = 150.0;       // max, min allowed angle, in degrees
var g_angleMin1  =  60.0;
                                //---------------
var g_angleNow2  =  0.0; 			  // Current rotation angle, in degrees.
var g_angleRate2 = -62.0;				// Rotation angle rate, in degrees/second.

                                //---------------
var g_posNow0 =  0.0;           // current position
var g_posRate0 = 0.6;           // position change rate, in distance/second.
var g_posMax0 =  0.5;           // max, min allowed for g_posNow;
var g_posMin0 = -0.5;           
                                // ------------------
var g_posNow1 =  0.0;           // current position
var g_posRate1 = 0.5;           // position change rate, in distance/second.
var g_posMax1 =  1.0;           // max, min allowed positions
var g_posMin1 = -1.0;
                                //---------------
                                //---------------
var g_angle01 = 0.0;                  // initial rotation angle
var g_angle01Rate = 45.0;           // rotation speed, in degrees/second 

var g_angle02 = 0.0;
var g_angle02Rate = 20.0;


var angleA = 0.0;
var angleB = 0.0;
var angleC = 0.0;


// For mouse/keyboard:------------------------
var g_show0 = 1;								//   "        "     Grid                              .
var g_show1 = 1;								// 	"					"			VBO1		"				"				" 
var g_show2 = 1;                //  "         "     Blin-phong    "       "       "


var g_isBlinn = false;
var g_fatt = 1;
var g_shiny;
var g_shinyInit;


var g_lightPosX;
var g_lightPosY;
var g_lightPosZ;

var g_lightDiffR;
var g_lightDiffG;
var g_lightDiffB;

var g_lightAmbiR;
var g_lightAmbiG;
var g_lightAmbiB;

var g_lightSpecR;
var g_lightSpecG;
var g_lightSpecB;

var g_lightPosXReset = document.getElementById('posX').value;
var g_lightPosYReset= document.getElementById('posY').value;
var g_lightPosZReset = document.getElementById('posZ').value;

var g_lightDiffRReset= document.getElementById('diffR').value;
var g_lightDiffGReset = document.getElementById('diffG').value;
var g_lightDiffBReset = document.getElementById('diffB').value;

var g_lightAmbiRReset = document.getElementById('ambiR').value;
var g_lightAmbiGReset = document.getElementById('ambiG').value;
var g_lightAmbiBReset = document.getElementById('ambiB').value;

var g_lightSpecRReset = document.getElementById('specR').value;
var g_lightSpecGReset = document.getElementById('specG').value;
var g_lightSpecBReset = document.getElementById('specB').value;

var g_EyeX = 6.5, g_EyeY = 5.5, g_EyeZ = 5.0;
var g_AtX = 1.0, g_AtY = 1.0, g_AtZ = 4.5;
var theta = 215;
var g_moveRate = 3.0;
var g_aimZDiff = g_AtZ - g_EyeZ;


// Key controls
window.addEventListener("keydown", myKeyDown, false);
window.addEventListener("keyup", myKeyUp, false);

function getUsrValues() {
  var usrPosX, usrPosY, usrPosZ, 
      usrAmbiR, usrAmbiG, usrAmbiB, 
      usrDiffR, usrDiffG, usrDiffB, 
      usrSpecR, usrSpecG, usrSpecB;

  // * Light position in world coords
  
  usrPosX = document.getElementById('posX').value;
  if(!isNaN(usrPosX)) g_lightPosX = usrPosX;

  usrPosY = document.getElementById('posY').value;
  if(!isNaN(usrPosY)) g_lightPosY = usrPosY

  usrPosZ = document.getElementById('posZ').value;
  if(!isNaN(usrPosZ)) g_lightPosZ = usrPosZ;

  // * Ambient light color

  usrAmbiR = document.getElementById('ambiR').value;
  if(!isNaN(usrAmbiR)) g_lightAmbiR = usrAmbiR;

  usrAmbiG = document.getElementById('ambiG').value;
  if(!isNaN(usrAmbiG)) g_lightAmbiG = usrAmbiG;

  usrAmbiB = document.getElementById('ambiB').value;
  if(!isNaN(usrAmbiB)) g_lightAmbiB = usrAmbiB;

  // * Diffuse light color

  usrDiffR = document.getElementById('diffR').value;
  if(!isNaN(usrDiffR)) g_lightDiffR = usrDiffR;

  usrDiffG = document.getElementById('diffG').value;
  if(!isNaN(usrDiffG)) g_lightDiffG = usrDiffG;

  usrDiffB = document.getElementById('diffB').value;
  if(!isNaN(usrDiffB)) g_lightDiffB = usrDiffB;

  // * Specular light color

  usrSpecR = document.getElementById('specR').value;
  if(!isNaN(usrSpecR)) g_lightSpecR = usrSpecR;

  usrSpecG = document.getElementById('specG').value;
  if(!isNaN(usrSpecG)) g_lightSpecG = usrSpecG;

  usrSpecB = document.getElementById('specB').value;
  if(!isNaN(usrSpecB)) g_lightSpecB = usrSpecB;
}


// GLOBAL CAMERA CONTROL:					// 
g_worldMat = new Matrix4();				// Changes CVV drawing axes to 'world' axes.
// (equivalently: transforms 'world' coord. numbers (x,y,z,w) to CVV coord. numbers)
// WHY?
// Lets mouse/keyboard functions set just one global matrix for 'view' and
// 'projection' transforms; then VBObox objects use it in their 'adjust()'
// member functions to ensure every VBObox draws its 3D parts and assemblies
// using the same 3D camera at the same 3D position in the same 3D world).


function main() {
//=============================================================================
  // Retrieve the HTML-5 <canvas> element where webGL will draw our pictures:
  g_canvasID = document.getElementById('webgl');	
  // Create the the WebGL rendering context: one giant JavaScript object that
  // contains the WebGL state machine adjusted by large sets of WebGL functions,
  // built-in variables & parameters, and member data. Every WebGL function call
  // will follow this format:  gl.WebGLfunctionName(args);

  // Create the the WebGL rendering context: one giant JavaScript object that
  // contains the WebGL state machine, adjusted by big sets of WebGL functions,
  // built-in variables & parameters, and member data. Every WebGL func. call
  // will follow this format:  gl.WebGLfunctionName(args);
  //SIMPLE VERSION:  gl = getWebGLContext(g_canvasID); 
  // Here's a BETTER version:
  gl = g_canvasID.getContext("webgl", { preserveDrawingBuffer: true});
	// This fancier-looking version disables HTML-5's default screen-clearing, so 
	// that our drawMain() 
	// function will over-write previous on-screen results until we call the 
	// gl.clear(COLOR_BUFFER_BIT); function. )
  if (!gl) {
    console.log('Failed to get the rendering context for WebGL');
    return;
  }
  gl.clearColor(0.2, 0.2, 0.2, 1);	  // RGBA color for clearing <canvas>

  gl.enable(gl.DEPTH_TEST);

  /*
//----------------SOLVE THE 'REVERSED DEPTH' PROBLEM:------------------------
  // IF the GPU doesn't transform our vertices by a 3D Camera Projection Matrix
  // (and it doesn't -- not until Project B) then the GPU will compute reversed 
  // depth values:  depth==0 for vertex z == -1;   (but depth = 0 means 'near') 
  //		    depth==1 for vertex z == +1.   (and depth = 1 means 'far').
  //
  // To correct the 'REVERSED DEPTH' problem, we could:
  //  a) reverse the sign of z before we render it (e.g. scale(1,1,-1); ugh.)
  //  b) reverse the usage of the depth-buffer's stored values, like this:
  gl.enable(gl.DEPTH_TEST); // enabled by default, but let's be SURE.

  gl.clearDepth(0.0);       // each time we 'clear' our depth buffer, set all
                            // pixel depths to 0.0  (1.0 is DEFAULT)
  gl.depthFunc(gl.GREATER); // draw a pixel only if its depth value is GREATER
                            // than the depth buffer's stored value.
                            // (gl.LESS is DEFAULT; reverse it!)
  //------------------end 'REVERSED DEPTH' fix---------------------------------
*/

  // Initialize each of our 'vboBox' objects: 
  worldBox.init(gl);		// VBO + shaders + uniforms + attribs for our 3D world,
                        // including ground-plane,                       
  gouraudBox.init(gl);		//  "		"		"  for 1st kind of shading & lighting
	phongBox.init(gl);    //  "   "   "  for 2nd kind of shading & lighting
	

  gl.clearColor(0.2, 0.2, 0.2, 1);	  // RGBA color for clearing <canvas>
  
  // ==============ANIMATION=============
  // Quick tutorials on synchronous, real-time animation in JavaScript/HTML-5: 
  //    https://webglfundamentals.org/webgl/lessons/webgl-animation.html
  //  or
  //  	http://creativejs.com/resources/requestanimationframe/
  //		--------------------------------------------------------
  // Why use 'requestAnimationFrame()' instead of the simpler-to-use
  //	fixed-time setInterval() or setTimeout() functions?  Because:
  //		1) it draws the next animation frame 'at the next opportunity' instead 
  //			of a fixed time interval. It allows your browser and operating system
  //			to manage its own processes, power, & computing loads, and to respond 
  //			to on-screen window placement (to skip battery-draining animation in 
  //			any window that was hidden behind others, or was scrolled off-screen)
  //		2) it helps your program avoid 'stuttering' or 'jittery' animation
  //			due to delayed or 'missed' frames.  Your program can read and respond 
  //			to the ACTUAL time interval between displayed frames instead of fixed
  //		 	fixed-time 'setInterval()' calls that may take longer than expected.
  //------------------------------------
  var tick = function() {		    // locally (within main() only), define our
    // self-calling animation function. 
    g_canvasID.width = innerWidth;
    g_canvasID.height = innerHeight * 0.70;
    
    
    requestAnimationFrame(tick, g_canvasID); // browser callback request; wait
    // til browser is ready to re-draw canvas, then
    timerAll();  // Update all time-varying params, and
    drawAll();                // Draw all the VBObox contents
    };
  //------------------------------------
  tick();                       // do it again!
}

function timerAll() {
//=============================================================================
// Find new values for all time-varying parameters used for on-screen drawing
  // use local variables to find the elapsed time.
  var nowMS = Date.now();             // current time (in milliseconds)
  var elapsedMS = nowMS - g_lastMS;   // 
  g_lastMS = nowMS;                   // update for next webGL drawing.
  if(elapsedMS > 1000.0) {            
    // Browsers won't re-draw 'canvas' element that isn't visible on-screen 
    // (user chose a different browser tab, etc.); when users make the browser
    // window visible again our resulting 'elapsedMS' value has gotten HUGE.
    // Instead of allowing a HUGE change in all our time-dependent parameters,
    // let's pretend that only a nominal 1/30th second passed:
    elapsedMS = 1000.0/30.0;
    }
  // Find new time-dependent parameters using the current or elapsed time:
  // Continuous rotation:
  g_angleNow0 = g_angleNow0 + (g_angleRate0 * elapsedMS) / 1000.0;
  g_angleNow1 = g_angleNow1 + (g_angleRate1 * elapsedMS) / 1000.0;
  g_angleNow2 = g_angleNow2 + (g_angleRate2 * elapsedMS) / 1000.0;
  g_angleNow0 %= 360.0;   // keep angle >=0.0 and <360.0 degrees  
  g_angleNow1 %= 360.0;   
  g_angleNow2 %= 360.0;
  // if(g_angleNow1 > g_angleMax1) { // above the max?
  //   g_angleNow1 = g_angleMax1;    // move back down to the max, and
  //   g_angleRate1 = -g_angleRate1; // reverse direction of change.
  //   }
  // else if(g_angleNow1 < g_angleMin1) {  // below the min?
  //   g_angleNow1 = g_angleMin1;    // move back up to the min, and
  //   g_angleRate1 = -g_angleRate1;
  //   }
  // Continuous movement:
  g_posNow0 += g_posRate0 * elapsedMS / 1000.0;
  g_posNow1 += g_posRate1 * elapsedMS / 1000.0;
  // apply position limits
  if(g_posNow0 > g_posMax0) {   // above the max?
    g_posNow0 = g_posMax0;      // move back down to the max, and
    g_posRate0 = -g_posRate0;   // reverse direction of change
    }
  else if(g_posNow0 < g_posMin0) {  // or below the min? 
    g_posNow0 = g_posMin0;      // move back up to the min, and
    g_posRate0 = -g_posRate0;   // reverse direction of change.
    }
  if(g_posNow1 > g_posMax1) {   // above the max?
    g_posNow1 = g_posMax1;      // move back down to the max, and
    g_posRate1 = -g_posRate1;   // reverse direction of change
    }
  else if(g_posNow1 < g_posMin1) {  // or below the min? 
    g_posNow1 = g_posMin1;      // move back up to the min, and
    g_posRate1 = -g_posRate1;   // reverse direction of change.
  }
  
var angleARate = 40.0;
var angleBRate = 20.0
var angleCRate = 20.0;

var angleAmax = 60.0;
var angleBmax = 50.0;
var angleCmax = 40.0;

  angleA = angleAmax * Math.sin(0.1 * Math.PI * angleARate * (performance.now()*0.0001));
  angleB = angleBmax * Math.sin(0.1 * Math.PI * angleBRate * (performance.now() * 0.0001));
  angleC = angleCmax * Math.sin(0.1 * Math.PI * angleCRate * (performance.now() * 0.0001));


}

function drawAll() {
//=============================================================================
  // Clear on-screen HTML-5 <canvas> object:
  // g_currMatl = getMatl(selectedInd1);
  matl1.setMatl(selectedInd1 + 1);
  matl2.setMatl(selectedInd2 + 1);
  g_shiny = document.getElementById("shiny").value;

  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  setCamera();
  getUsrValues();
  // console.log("angles: ", angleA, angleB, angleC)
  
var b4Draw = Date.now();
var b4Wait = b4Draw - g_lastMS;
  // set material

	if(g_show0 == 1) {	// IF user didn't press HTML button to 'hide' VBO0:
	  worldBox.switchToMe();  // Set WebGL to render from this VBObox.
		worldBox.adjust();		  // Send new values for uniforms to the GPU, and
		worldBox.draw();			  // draw our VBO's contents using our shaders.
  }

  if(g_show1 == 1) { // IF user didn't press HTML button to 'hide' VBO1:
    gouraudBox.switchToMe();  // Set WebGL to render from this VBObox.
  	gouraudBox.adjust();		  // Send new values for uniforms to the GPU, and
    gouraudBox.draw();			  // draw our VBO's contents using our shaders.
  }
  else if (g_show1 == 0) {  // IF user didn't press HTML button to 'hide' VBO1:
    phongBox.switchToMe();  // Set WebGL to render from this VBObox.
  	phongBox.adjust();		  // Send new values for uniforms to the GPU, and
  	phongBox.draw();			  // draw our VBO's contents using our shaders.
  }

	// if(g_show2 == 1) { // IF user didn't press HTML button to 'hide' VBO2:
	//   phongBox.switchToMe();  // Set WebGL to render from this VBObox.
  // 	phongBox.adjust();		  // Send new values for uniforms to the GPU, and
  // 	phongBox.draw();			  // draw our VBO's contents using our shaders.
  // 	}
/* // ?How slow is our own code?  	
var aftrDraw = Date.now();
var drawWait = aftrDraw - b4Draw;
console.log("wait b4 draw: ", b4Wait, "drawWait: ", drawWait, "mSec");
*/
}

function VBO0toggle() {
//=============================================================================
// Called when user presses HTML-5 button 'Show/Hide VBO0'.
  if(g_show0 != 1) g_show0 = 1;				// show,
  else g_show0 = 0;										// hide.
  console.log('g_show0: '+g_show0);
}

function GPtoggle() {
//=============================================================================
// Called when user presses HTML-5 button
  if (g_show1 != 1) {
    g_show1 = 1;			// show,
    document.getElementById('shading').innerHTML = "Phong";
  }
  else {
    g_show1 = 0;
    document.getElementById('shading').innerHTML = "Gouraud";
  }// hide.
  console.log('g_show1: ' + g_show1);
}

function toggleBlinn() {
  if (!g_isBlinn) {
    g_isBlinn = true;
    document.getElementById('lighting').innerHTML = "Blinn";
  }
  else {
    g_isBlinn = false;
    document.getElementById('lighting').innerHTML = "Phong";
  }
  console.log('g_isBlinn: '+ g_isBlinn);
}


function toRadians(angle) {
	return angle * (Math.PI/180);
}
function setCamera() {
//============================================================================
// PLACEHOLDER:  sets a fixed camera at a fixed position for use by
// ALL VBObox objects.  REPLACE This with your own camera-control code.

  g_worldMat.setIdentity();
  gl.viewport(0,											 				// Viewport lower-left corner
  0, 			// location(in pixels)
  g_canvasID.width, 					// viewport width,
  g_canvasID.height);			// viewport height in pixels.


	g_worldMat.perspective(40.0,   // FOVY: top-to-bottom vertical image angle, in degrees
  g_canvasID.width/g_canvasID.height,   // Image Aspect Ratio: camera lens width/height
                      1.0,   // camera z-near distance (always positive; frustum begins at z = -znear)
                      100.0);  // camera z-far distance (always positive; frustum ends at z = -zfar)
  
  g_AtX = g_EyeX + Math.cos(toRadians(theta));
  g_AtY = g_EyeY + Math.sin(toRadians(theta));
  g_AtZ = g_EyeZ + g_aimZDiff;
                    
  g_worldMat.lookAt(g_EyeX, g_EyeY, g_EyeZ, // eye position
                    g_AtX, g_AtY, g_AtZ,                  // look-at point 
                    0, 0, 1);   
	
  // READY to draw in the 'world' coordinate system.
//------------END COPY

}

function myKeyDown(ev) {
	var xd = g_EyeX - g_AtX;
	var yd = g_EyeY - g_AtY;
	var zd = g_EyeZ - g_AtZ;

	var lxy = Math.sqrt(xd * xd + yd * yd);
	var l = Math.sqrt(xd * xd + yd * yd + zd * zd);

  var moveRateRad = toRadians(g_moveRate);

	switch (ev.keyCode) { // keyboard listener
		// ---- move camera around
		case 87:    // w
			g_AtX = g_AtX - 0.1 * (xd / l);
			g_AtY = g_AtY - 0.1 * (yd / l);
			g_AtZ = g_AtZ - 0.1 * (zd / l);

			g_EyeX = g_EyeX - 0.1 * (xd / l);
			g_EyeY = g_EyeY - 0.1 * (yd / l);
			g_EyeZ = g_EyeZ - 0.1 * (zd / l);
			break;

		case 83:    // s
			g_AtX = g_AtX + 0.1 * (xd / l);
			g_AtY = g_AtY + 0.1 * (yd / l);
			g_AtZ = g_AtZ + 0.1 * (zd / l);

			g_EyeX = g_EyeX + 0.1 * (xd / l);
			g_EyeY = g_EyeY + 0.1 * (yd / l);
			g_EyeZ = g_EyeZ + 0.1 * (zd / l);
			break;

		case 68:    // a
			g_EyeX = g_EyeX - 0.1 * yd / lxy;
			g_EyeY = g_EyeY + 0.1 * xd / lxy;
			g_AtX -= 0.1 * yd / lxy;
			g_AtY += 0.1 * xd / lxy;

			break;
		case 65:    // d
			g_EyeX = g_EyeX + 0.1 * yd / lxy;
			g_EyeY = g_EyeY - 0.1 * xd / lxy;
			g_AtX += 0.1 * yd / lxy;
			g_AtY -= 0.1 * xd / lxy;
			break;
		
		// tilt camera
		case 74:   // j
			theta += g_moveRate;

			if(theta > 360) theta -= 360.0;
			if(theta < 0) theta += 360.0;
			break;
		case 76:   //l
      theta -= g_moveRate;
			if(theta > 360) theta -= 360.0;
      if (theta < 0) theta += 360.0;
			break;
    case 73: //i
        g_aimZDiff += moveRateRad;
        break;
    case 75:  // k
        g_aimZDiff -= moveRateRad;
        break;

		default:
        	console.log('myKeyDown()--keycode=', ev.keyCode, ', charCode=', ev.charCode);
     break;
    }

}

function myKeyUp(ev) {
	//===============================================================================
	// Called when user releases ANY key on the keyboard; captures scancodes well

	console.log('myKeyUp()--keyCode=' + ev.keyCode + ' released.');
}
  
function getMatl(selectedIndex) {
  if ( g_currMatl != selectedIndex+1){
    var matl = new Material(selectedIndex+1);
    document.getElementById('shiny').value = g_shinyInit;
    console.log("new chosen", selectedIndex)
    return selectedIndex + 1;
  }
  return g_currMatl;
}

function resetShiny() {
  document.getElementById('shiny').value = g_shinyInit;
}

function att (num) {
  g_fatt = num;
  }