//3456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_
// (JT: why the numbers? counts columns, helps me keep 80-char-wide listings)

// Tabs set to 2

/*=====================
  VBObox-Lib.js library: 
  ===================== 
Note that you don't really need 'VBObox' objects for any simple, 
    beginner-level WebGL/OpenGL programs: if all vertices contain exactly 
		the same attributes (e.g. position, color, surface normal), and use 
		the same shader program (e.g. same Vertex Shader and Fragment Shader), 
		then our textbook's simple 'example code' will suffice.
		  
***BUT*** that's rare -- most genuinely useful WebGL/OpenGL programs need 
		different sets of vertices with  different sets of attributes rendered 
		by different shader programs.  THUS a customized VBObox object for each 
		VBO/shader-program pair will help you remember and correctly implement ALL 
		the WebGL/GLSL steps required for a working multi-shader, multi-VBO program.
		
One 'VBObox' object contains all we need for WebGL/OpenGL to render on-screen a 
		set of shapes made from vertices stored in one Vertex Buffer Object (VBO), 
		as drawn by calls to one 'shader program' that runs on your computer's 
		Graphical Processing Unit(GPU), along with changes to values of that shader 
		program's one set of 'uniform' varibles.  
The 'shader program' consists of a Vertex Shader and a Fragment Shader written 
		in GLSL, compiled and linked and ready to execute as a Single-Instruction, 
		Multiple-Data (SIMD) parallel program executed simultaneously by multiple 
		'shader units' on the GPU.  The GPU runs one 'instance' of the Vertex 
		Shader for each vertex in every shape, and one 'instance' of the Fragment 
		Shader for every on-screen pixel covered by any part of any drawing 
		primitive defined by those vertices.
The 'VBO' consists of a 'buffer object' (a memory block reserved in the GPU),
		accessed by the shader program through its 'attribute' variables. Shader's
		'uniform' variable values also get retrieved from GPU memory, but their 
		values can't be changed while the shader program runs.  
		Each VBObox object stores its own 'uniform' values as vars in JavaScript; 
		its 'adjust()'	function computes newly-updated values for these uniform 
		vars and then transfers them to the GPU memory for use by shader program.
EVENTUALLY you should replace 'cuon-matrix-quat03.js' with the free, open-source
   'glmatrix.js' library for vectors, matrices & quaternions: Google it!
		This vector/matrix library is more complete, more widely-used, and runs
		faster than our textbook's 'cuon-matrix-quat03.js' library.  
		--------------------------------------------------------------
		I recommend you use glMatrix.js instead of cuon-matrix-quat03.js
		--------------------------------------------------------------
		for all future WebGL programs. 
You can CONVERT existing cuon-matrix-based programs to glmatrix.js in a very 
    gradual, sensible, testable way:
		--add the glmatrix.js library to an existing cuon-matrix-based program;
			(but don't call any of its functions yet).
		--comment out the glmatrix.js parts (if any) that cause conflicts or in	
			any way disrupt the operation of your program.
		--make just one small local change in your program; find a small, simple,
			easy-to-test portion of your program where you can replace a 
			cuon-matrix object or function call with a glmatrix function call.
			Test; make sure it works. Don't make too large a change: it's hard to fix!
		--Save a copy of this new program as your latest numbered version. Repeat
			the previous step: go on to the next small local change in your program
			and make another replacement of cuon-matrix use with glmatrix use. 
			Test it; make sure it works; save this as your next numbered version.
		--Continue this process until your program no longer uses any cuon-matrix
			library features at all, and no part of glmatrix is commented out.
			Remove cuon-matrix from your library, and now use only glmatrix.

	------------------------------------------------------------------
	VBObox -- A MESSY SET OF CUSTOMIZED OBJECTS--NOT REALLY A 'CLASS'
	------------------------------------------------------------------
As each 'VBObox' object can contain:
  -- a DIFFERENT GLSL shader program, 
  -- a DIFFERENT set of attributes that define a vertex for that shader program, 
  -- a DIFFERENT number of vertices to used to fill the VBOs in GPU memory, and 
  -- a DIFFERENT set of uniforms transferred to GPU memory for shader use.  
  THUS:
		I don't see any easy way to use the exact same object constructors and 
		prototypes for all VBObox objects.  Every additional VBObox objects may vary 
		substantially, so I recommend that you copy and re-name an existing VBObox 
		prototype object, and modify as needed, as shown here. 
		(e.g. to make the VBObox3 object, copy the VBObox2 constructor and 
		all its prototype functions, then modify their contents for VBObox3 
		activities.)

*/

// Written for EECS 351-2,	Intermediate Computer Graphics,
//							Northwestern Univ. EECS Dept., Jack Tumblin
// 2016.05.26 J. Tumblin-- Created; tested on 'TwoVBOs.html' starter code.
// 2017.02.20 J. Tumblin-- updated for EECS 351-1 use for Project C.
// 2018.04.11 J. Tumblin-- minor corrections/renaming for particle systems.
//    --11e: global 'gl' replaced redundant 'myGL' fcn args;
//    --12: added 'SwitchToMe()' fcn to simplify 'init()' function and to fix
//      weird subtle errors that sometimes appear when we alternate 'adjust()'
//      and 'draw()' functions of different VBObox objects. CAUSE: found that
//      only the 'draw()' function (and not the 'adjust()' function) made a full
//      changeover from one VBObox to another; thus calls to 'adjust()' for one
//      VBObox could corrupt GPU contents for another.
//      --Created vboStride, vboOffset members to centralize VBO layout in the
//      constructor function.
//    -- 13 (abandoned) tried to make a 'core' or 'resuable' VBObox object to
//      which we would add on new properties for shaders, uniforms, etc., but
//      I decided there was too little 'common' code that wasn't customized.
//=============================================================================


//=============================================================================
//=============================================================================
var floatsPerVertex = 6;

function makeGroundGrid() {
	//==============================================================================
	// Create a list of vertices that create a large grid of lines in the x,y plane
	// centered at the origin.  Draw this shape using the GL_LINES primitive.
	
		var xcount = 100;			// # of lines to draw in x,y to make the grid.
		var ycount = 100;		
    var xymax = 20.0;
  
  // grid size; extends to cover +/-xymax in x and y.
  var xColr = new Float32Array([1.0, 1.0, 0.3]);	// bright yellow
	var yColr = new Float32Array([0.5, 1.0, 0.5]);	// bright green.
		
	// Create an (global) array to hold this ground-plane's vertices:
	gndVerts = new Float32Array(floatsPerVertex*2*(xcount+ycount));
							// draw a grid made of xcount+ycount lines; 2 vertices per line.
							
	var xgap = xymax/(xcount-1);		// HALF-spacing between lines in x,y;
	var ygap = xymax/(ycount-1);		// (why half? because v==(0line number/2))
		
		// First, step thru x values as we make vertical lines of constant-x:
		for(v=0, j=0; v<2*xcount; v++, j+= floatsPerVertex) {
			if(v%2==0) {	// put even-numbered vertices at (xnow, -xymax, 0)
				gndVerts[j  ] = -xymax + (v  )*xgap;	// x
				gndVerts[j+1] = -xymax;								// y
				gndVerts[j+2] = 0.0;									// z
			}
			else {				// put odd-numbered vertices at (xnow, +xymax, 0).
				gndVerts[j  ] = -xymax + (v-1)*xgap;	// x
				gndVerts[j+1] = xymax;								// y
				gndVerts[j+2] = 0.0;									// z
			}
			gndVerts[j+3] = xColr[0];			// red
			gndVerts[j+4] = xColr[1];			// grn
			gndVerts[j+5] = xColr[2];			// blu
		}
		// Second, step thru y values as wqe make horizontal lines of constant-y:
		// (don't re-initialize j--we're adding more vertices to the array)
		for(v=0; v<2*ycount; v++, j+= floatsPerVertex) {
			if(v%2==0) {		// put even-numbered vertices at (-xymax, ynow, 0)
				gndVerts[j  ] = -xymax;								// x
				gndVerts[j+1] = -xymax + (v  )*ygap;	// y
				gndVerts[j+2] = 0.0;									// z
			}
			else {					// put odd-numbered vertices at (+xymax, ynow, 0).
				gndVerts[j  ] = xymax;								// x
				gndVerts[j+1] = -xymax + (v-1)*ygap;	// y
				gndVerts[j+2] = 0.0;									// z
			}
			gndVerts[j+3] = yColr[0];			// red
			gndVerts[j+4] = yColr[1];			// grn
			gndVerts[j+5] = yColr[2];			// blu
		}
}
function makeSphere() {
  //==============================================================================
  // Make a sphere from one OpenGL TRIANGLE_STRIP primitive.   Make ring-like 
  // equal-lattitude 'slices' of the sphere (bounded by planes of constant z), 
  // and connect them as a 'stepped spiral' design (see makeCylinder) to build the
  // sphere from one triangle strip.
  var slices = 41;		// # of slices of the sphere along the z axis. >=3 req'd
  // (choose odd # or prime# to avoid accidental symmetry)
  var sliceVerts = 41;	// # of vertices around the top edge of the slice
  // (same number of vertices on bottom of slice, too)
  var topColr = new Float32Array([0.5, 0.5, 0.5]);	// North Pole:
  var equColr = new Float32Array([.3, .3, .3]);	// Equator:    
  var botColr = new Float32Array([1, 1, 1]);	// South Pole: 
  var sliceAngle = Math.PI / slices;	// lattitude angle spanned by one slice.

  // Create a (global) array to hold this sphere's vertices:
  sphVerts = new Float32Array(((slices * 2 * sliceVerts) - 2) * floatsPerVertex);
  // # of vertices * # of elements needed to store them. 
  // each slice requires 2*sliceVerts vertices except 1st and
  // last ones, which require only 2*sliceVerts-1.

  // Create dome-shaped top slice of sphere at z=+1
  // s counts slices; v counts vertices; 
  // j counts array elements (vertices * elements per vertex)
  var cos0 = 0.0;					// sines,cosines of slice's top, bottom edge.
  var sin0 = 0.0;
  var cos1 = 0.0;
  var sin1 = 0.0;
  var j = 0;							// initialize our array index
  var isLast = 0;
  var isFirst = 1;
  for (s = 0; s < slices; s++) {	// for each slice of the sphere,
      // find sines & cosines for top and bottom of this slice
      if (s == 0) {
          isFirst = 1;	// skip 1st vertex of 1st slice.
          cos0 = 1.0; 	// initialize: start at north pole.
          sin0 = 0.0;
      }
      else {					// otherwise, new top edge == old bottom edge
          isFirst = 0;
          cos0 = cos1;
          sin0 = sin1;
      }								// & compute sine,cosine for new bottom edge.
      cos1 = Math.cos((s + 1) * sliceAngle);
      sin1 = Math.sin((s + 1) * sliceAngle);
      // go around the entire slice, generating TRIANGLE_STRIP verts
      // (Note we don't initialize j; grows with each new attrib,vertex, and slice)
      if (s == slices - 1) isLast = 1;	// skip last vertex of last slice.
      for (v = isFirst; v < 2 * sliceVerts - isLast; v++, j += floatsPerVertex) {
          if (v % 2 == 0) {				// put even# vertices at the the slice's top edge
              // (why PI and not 2*PI? because 0 <= v < 2*sliceVerts
              // and thus we can simplify cos(2*PI(v/2*sliceVerts))  
              sphVerts[j] = sin0 * Math.cos(Math.PI * (v) / sliceVerts);
              sphVerts[j + 1] = sin0 * Math.sin(Math.PI * (v) / sliceVerts);
              sphVerts[j + 2] = cos0;
          }
          else { 	// put odd# vertices around the slice's lower edge;
              // x,y,z,w == cos(theta),sin(theta), 1.0, 1.0
              // 					theta = 2*PI*((v-1)/2)/capVerts = PI*(v-1)/capVerts
              sphVerts[j] = sin1 * Math.cos(Math.PI * (v - 1) / sliceVerts);		// x
              sphVerts[j + 1] = sin1 * Math.sin(Math.PI * (v - 1) / sliceVerts);		// y
              sphVerts[j + 2] = cos1;			
          }
          if (s == 0) {	// finally, set some interesting colors for vertices:

              sphVerts[j + 3] = sin1 * Math.cos(Math.PI * (v - 1) / sliceVerts);
              sphVerts[j + 4] = sin1 * Math.sin(Math.PI * (v - 1) / sliceVerts);
              sphVerts[j + 5] = cos1;
          }
          else if (s == slices - 1) {
              sphVerts[j + 3] = sin1 * Math.cos(Math.PI * (v - 1) / sliceVerts);
              sphVerts[j + 4] = sin1 * Math.sin(Math.PI * (v - 1) / sliceVerts);
              sphVerts[j + 5] = cos1;
          }
          else {
              sphVerts[j + 3] = sin1 * Math.cos(Math.PI * (v - 1) / sliceVerts);
              sphVerts[j + 4] = sin1 * Math.sin(Math.PI * (v - 1) / sliceVerts);
              sphVerts[j + 5] = cos1;
          }

      }
  }
}
function makeCube() {

  cubeVerts = new Float32Array([
  // Vertex coordinates(x,y,z,w) and color (R,G,B) for a color tetrahedron:
  //    Apex on +z axis; equilateral triangle base at z=0

    // +x face: RED
     1.0, -0.0, -0.2,  1, 0, 0,// Node 3
     1.0, 2.0, -0.2, 1, 0, 0,// Node 2
     1.0, 2.0, 0.2, 1, 0, 0,// Node 4

     1.0, 2.0, 0.2, 1, 0, 0,// Node 4
     1.0, -0.0, 0.2, 1, 0, 0,// Node 7
     1.0, -0.0, -0.2, 1, 0, 0,// Node 3

    // +y face: GREEN
    -1.0, 2.0, -0.2, 0, 1, 0,// Node 1
    -1.0, 2.0, 0.2, 0, 1, 0,// Node 5
     1.0, 2.0, 0.2, 0, 1, 0,// Node 4

     1.0, 2.0, 0.2, 0, 1, 0,// Node 4
     1.0, 2.0, -0.2, 0, 1, 0,// Node 2 
    -1.0, 2.0, -0.2, 0, 1, 0,// Node 1

    // +z face: BLUE
    -1.0, 2.0, 0.2, 0, 0, 1,// Node 5
    -1.0, -0.0, 0.2, 0, 0, 1,// Node 6
     1.0, -0.0, 0.2, 0, 0, 1,// Node 7

     1.0, -0.0, 0.2, 0, 0, 1,// Node 7
     1.0, 2.0, 0.2,  0, 0, 1,// Node 4
    -1.0, 2.0, 0.2,  0, 0, 1,// Node 5

    // -x face: CYAN
    -1.0, -0.0, 0.2,-1, 0, 0,// Node 6 
    -1.0, 2.0, 0.2,  -1, 0, 0,// Node 5 
    -1.0, 2.0, -0.2, -1, 0, 0,// Node 1

    -1.0, 2.0, -0.2, -1, 0, 0,// Node 1
    -1.0, -0.0, -0.2, -1, 0, 0,// Node 0  
    -1.0, -0.0, 0.2, -1, 0, 0,// Node 6  

    // -y face: MAGENTA
     1.0, -0.0, -0.2, 0, -1, 0,// Node 3
     1.0, -0.0, 0.2, 0, -1, 0,// Node 7
    -1.0, -0.0, 0.2, 0, -1, 0,// Node 6

    -1.0, -0.0, 0.2, 0, -1, 0,// Node 6
    -1.0, -0.0, -0.2, 0, -1, 0,// Node 0
     1.0, -0.0, -0.2, 0, -1, 0,// Node 3

     // -z face: YELLOW
     1.0, 2.0, -0.2, 0, 0, -1,// Node 2
     1.0, -0.0, -0.2, 0, 0, -1,// Node 3
    -1.0, -0.0, -0.2,0, 0, -1,// Node 0   

    -1.0, -0.0, -0.2, 0, 0, -1,// Node 0
    -1.0, 2.0, -0.2,  0, 0, -1,// Node 1
     1.0, 2.0, -0.2,  0, 0, -1,// Node 2

  ]);

}
function makeCylinder(topRadius, bottomRadius) {
  //==============================================================================
  // Make a cylinder shape from one TRIANGLE_STRIP drawing primitive, using the
  // 'stepped spiral' design (Method 2) described in the class lecture notes.
  // Cylinder center at origin, encircles z axis, radius 1, top/bottom at z= +/-1.
  //
  
  var topColr = new Float32Array([0.8, 0.8, 0.0]);	// light yellow top,
  var walColr = new Float32Array([0.2, 0.6, 0.2]);	// dark green walls,
  var botColr = new Float32Array([0.2, 0.3, 0.7]);	// light blue bottom,
  var ctrColr = new Float32Array([0.1, 0.1, 0.1]); // near black end-cap centers,
  var errColr = new Float32Array([1.0, 0.2, 0.2]);	// Bright-red trouble color.
  
  var capVerts = 6;	// # of vertices around the topmost 'cap' of the shape
  //var topRadius = 0.9;		// radius of top of cylinder (bottom is always 1.0)
  
  // Create a (global) array to hold all of this cylinder's vertices;
  cylVerts = new Float32Array(  ((capVerts*6) -2) * floatsPerVertex);
  // # of vertices * # of elements needed to store them. How many vertices?
                      // Cylinder bottom end-cap:   (2*capVerts) -1  vertices;
                      // (includes blue transition-edge that links end-cap & wall);
                      // + Cylinder wall requires   (2*capVerts) vertices;
                      // + Cylinder top end-cap:    (2*capVerts) -1 vertices
                      // (includes green transition-edge that links wall & endcap).
  
    // Create circle-shaped bottom cap of cylinder at z=-1.0, radius 1.0,
    // with (capVerts*2)-1 vertices, BUT STOP before you create it's last vertex.
    // That last vertex forms the 'transition' edge from the bottom cap to the 
    // wall (shown in blue in lecture notes), & we make it in the next for() loop.
    // 
    // v counts vertices: j counts array elements (vertices * elements per vertex)
    for(v=0,j=0;   v<(2*capVerts)-1;   v++,j+=floatsPerVertex) {	
      // START at vertex v = 0; on x-axis on end-cap's outer edge, at xyz = 1,0,-1.
      // END at the vertex 2*(capVerts-1), the last outer-edge vertex before 
      // we reach the starting vertex at 1,0,-1. 
      if(v%2 ==0)
      {				// put even# vertices around bottom cap's outer edge,starting at v=0.
              // visit each outer-edge location only once; don't return to 
              // to the location of the v=0 vertex (at 1,0,-1).
              // x,y,z,w == cos(theta),sin(theta),-1.0, 1.0, 
              // 		where	theta = 2*PI*((v/2)/capVerts) = PI*v/capVerts
        cylVerts[j  ] =bottomRadius* Math.cos(Math.PI*v/capVerts);			// x
        cylVerts[j+1] = -0.4;			// y
        //	(Why not 2*PI? because 0 < =v < 2*capVerts,
        //	 so we can simplify cos(2*PI * (v/(2*capVerts))
        cylVerts[j+2] =bottomRadius*Math.sin(Math.PI*v/capVerts);;	// z
        // r,g,b = botColr[] 
        cylVerts[j+3]=botColr[0]; 
        cylVerts[j+4]=botColr[1]; 
        cylVerts[j+5]=botColr[2];
      }
      else {	// put odd# vertices at center of cylinder's bottom cap:
        cylVerts[j  ] = 0.0; 			// x,y,z,w == 0,0,-1,1; centered on z axis at -1.
        cylVerts[j+1] = -0.4;	
        cylVerts[j+2] =0.0; 
        cylVerts[j+3]=ctrColr[0]; 
        cylVerts[j+4]=ctrColr[1]; 
        cylVerts[j+5]=ctrColr[2];
      }
    }
    // Create the cylinder side walls, made of 2*capVerts vertices.
    // v counts vertices within the wall; j continues to count array elements
    // START with the vertex at 1,0,-1 (completes the cylinder's bottom cap;
    // completes the 'transition edge' drawn in blue in lecture notes).
    for(v=0; v< 2*capVerts;   v++, j+=floatsPerVertex) {
      if(v%2==0)	// count verts from zero again, 
                  // and put all even# verts along outer edge of bottom cap:
      {		
          cylVerts[j  ] = bottomRadius*Math.cos(Math.PI*(v)/capVerts);		// x
          cylVerts[j+1] = -0.4;		// y
          cylVerts[j+2] =bottomRadius*Math.sin(Math.PI*(v)/capVerts);	// ==z  BOTTOM cap,
          // r,g,b = walColr[]				
          cylVerts[j+3]=walColr[0]; 
          cylVerts[j+4]=walColr[1]; 
          cylVerts[j+5]=walColr[2];			
        if(v==0) {		// UGLY TROUBLESOME vertex--shares its 1 color with THREE
                      // triangles; 1 in cap, 1 in step, 1 in wall.
            cylVerts[j+3] = errColr[0]; 
            cylVerts[j+4] = errColr[1];
            cylVerts[j+5] = errColr[2];		// (make it red; see lecture notes)
          }
      }
      else		// position all odd# vertices along the top cap (not yet created)
      {
          cylVerts[j  ] = topRadius * Math.cos(Math.PI*(v-1)/capVerts);		// x
          cylVerts[j+1] = 0.4;		// y
          cylVerts[j+2] = topRadius * Math.sin(Math.PI*(v-1)/capVerts);	// == z TOP cap,

          // r,g,b = walColr;
          cylVerts[j+3]=walColr[0]; 
          cylVerts[j+4]=walColr[1]; 
          cylVerts[j+5]=walColr[2];			
      }
    }
    // Complete the cylinder with its top cap, made of 2*capVerts -1 vertices.
    // v counts the vertices in the cap; j continues to count array elements.
    for(v=0; v < (2*capVerts -1); v++, j+= floatsPerVertex) {
      // count vertices from zero again, and
      if(v%2==0) {	// position even #'d vertices around top cap's outer edge.
        cylVerts[j  ] = topRadius * Math.cos(Math.PI*(v)/capVerts);		// x
        cylVerts[j+1] = 0.4;		// y
        cylVerts[j + 2] = topRadius * Math.sin(Math.PI * (v) / capVerts);	// z
      
        // r,g,b = topColr[]
        cylVerts[j+3]=topColr[0]; 
        cylVerts[j+4]=topColr[1]; 
        cylVerts[j+5]=topColr[2];
        if(v==0) {	// UGLY TROUBLESOME vertex--shares its 1 color with THREE
                      // triangles; 1 in cap, 1 in step, 1 in wall.
            cylVerts[j+3] = errColr[0]; 
            cylVerts[j+4] = errColr[1];
            cylVerts[j+5] = errColr[2];		// (make it red; see lecture notes)
        }		
      }
      else {				// position odd#'d vertices at center of the top cap:
        cylVerts[j  ] = 0.0; 			// x,y,z,w == 0,0,-1,1
        cylVerts[j+1] = 0.4;	
        cylVerts[j + 2] = 0.0; 
        
        // r,g,b = topColr[]
        cylVerts[j+3]=ctrColr[0]; 
        cylVerts[j+4]=ctrColr[1]; 
        cylVerts[j+5]=ctrColr[2];
      }
  }
}

function makeAxes() {
	axesVerts =new Float32Array([
		// X-axis
		-30.0, 0.0, 0.0,     1.0,  0.0,  0.0,  // v0 Red
		30.0, 0.0, 0.0,     1.0,  0.0,  0.0,  // v1 Red

		// Y-axis
		0.0, -30.0, 0.0,    0.0,  0.0,  1.0,  // v0 Blue
		0.0, 30.0, 0.0,     0.0,  0.0,  1.0,  // v1 Blue

		// Z-Axis
		0.0, 0.0, -30.0,    0.0,  1.0,  0.0,  // v0 Green
		0.0, 0.0, 30.0,     0.0,  1.0,  0.0,  // v1 Green
	]);
	
}


function VBObox0() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox0' object that holds all data and fcns
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate shader program (a vertex-shader & fragment-shader pair) and one
// set of 'uniform' variables.

// Constructor goal: 
// Create and set member vars that will ELIMINATE ALL LITERALS (numerical values 
// written into code) in all other VBObox functions. Keeping all these (initial)
// values here, in this one coonstrutor function, ensures we can change them 
// easily WITHOUT disrupting any other code, ever!
  
	this.VERT_SRC =	//--------------------- VERTEX SHADER source code 
 `precision highp float;				// req'd in OpenGL ES if we use 'float'
  //
  uniform mat4 u_ModelMat0;
  attribute vec4 a_Pos0;
  attribute vec3 a_Colr0;
  varying vec3 v_Colr0;
  //
  void main() {
    gl_Position = u_ModelMat0 * a_Pos0;
  	 v_Colr0 = a_Colr0;
   }`;

	this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code 
 `precision mediump float;
  varying vec3 v_Colr0;
  void main() {
    gl_FragColor = vec4(v_Colr0, 1.0);
  }`;


  makeGroundGrid();
  makeAxes();

  this.vboSize = (axesVerts.length + gndVerts.length);	
  // this.vboSize = (axesVerts.length);	
  this.vboVerts = this.vboSize / floatsPerVertex;	
  
	this.vboContents = //---------------------------------------------------------
            new Float32Array(this.vboSize);			// Array of vertex attribute values we will
  															        // transfer to GPU's vertex buffer object (VBO)
	
  // Copy them:  remember where to start for each shape:
	GroundStart = 0;
	for(i=0,j=0; j< gndVerts.length; i++,j++) {
      this.vboContents[i] = gndVerts[j];
  }
  
	axesStart = i;
	for(j=0; j< axesVerts.length; i++, j++) {
		this.vboContents[i] = axesVerts[j];
  }	
  
	// axesStart = 0;
	// for(i=0, j=0; j< axesVerts.length; i++, j++) {
	// 	this.vboContents[i] = axesVerts[j];
	// }	
  console.log(this.vboContents.length)



					// # of vertices held in 'vboContents' array
	this.FSIZE = this.vboContents.BYTES_PER_ELEMENT;
	                              // bytes req'd by 1 vboContents array element;
																// (why? used to compute stride and offset 
																// in bytes for vertexAttribPointer() calls)
  this.vboBytes = this.vboContents.length * this.FSIZE;               
                                // total number of bytes stored in vboContents
                                // (#  of floats in vboContents array) * 
                                // (# of bytes/float).
  console.log(this.vboBytes, this.vboVerts)
  
	this.vboStride = this.vboBytes / this.vboVerts; 
	                              // (== # of bytes to store one complete vertex).
	                              // From any attrib in a given vertex in the VBO, 
	                              // move forward by 'vboStride' bytes to arrive 
	                              // at the same attrib for the next vertex. 
                                

	            //----------------------Attribute sizes
  this.vboFcount_a_Pos0 =  3;    // # of floats in the VBO needed to store the
                                // attribute named a_Pos0. (4: x,y,z values)
  this.vboFcount_a_Colr0 = 3;   // # of floats for this attrib (r,g,b values)

  console.assert((this.vboFcount_a_Pos0 +     // check the size of each and
                  this.vboFcount_a_Colr0) *   // every attribute in our VBO
                  this.FSIZE == this.vboStride, // for agreeement with'stride'
                  "Uh oh! VBObox0.vboStride disagrees with attribute-size values!");

              //----------------------Attribute offsets  
	this.vboOffset_a_Pos0 = 0;    // # of bytes from START of vbo to the START
	                              // of 1st a_Pos0 attrib value in vboContents[]
  this.vboOffset_a_Colr0 = this.vboFcount_a_Pos0 * this.FSIZE;    
                                // (4 floats * bytes/float) 
                                // # of bytes from START of vbo to the START
                                // of 1st a_Colr0 attrib value in vboContents[]
	            //-----------------------GPU memory locations:
	this.vboLoc;									// GPU Location for Vertex Buffer Object, 
	                              // returned by gl.createBuffer() function call
	this.shaderLoc;								// GPU Location for compiled Shader-program  
	                            	// set by compile/link of VERT_SRC and FRAG_SRC.
								          //------Attribute locations in our shaders:
	this.a_PosLoc;								// GPU location for 'a_Pos0' attribute
	this.a_ColrLoc;								// GPU location for 'a_Colr0' attribute

	            //---------------------- Uniform locations &values in our shaders
	this.ModelMat = new Matrix4();	// Transforms CVV axes to model axes.
	this.u_ModelMatLoc;							// GPU location for u_ModelMat uniform
}

VBObox0.prototype.init = function() {
//=============================================================================
// Prepare the GPU to use all vertices, GLSL shaders, attributes, & uniforms 
// kept in this VBObox. (This function usually called only once, within main()).
// Specifically:
// a) Create, compile, link our GLSL vertex- and fragment-shaders to form an 
//  executable 'program' stored and ready to use inside the GPU.  
// b) create a new VBO object in GPU memory and fill it by transferring in all
//  the vertex data held in our Float32array member 'VBOcontents'. 
// c) Find & save the GPU location of all our shaders' attribute-variables and 
//  uniform-variables (needed by switchToMe(), adjust(), draw(), reload(), etc.)
// -------------------
// CAREFUL!  before you can draw pictures using this VBObox contents, 
//  you must call this VBObox object's switchToMe() function too!
//--------------------
// a) Compile,link,upload shaders-----------------------------------------------
	this.shaderLoc = createProgram(gl, this.VERT_SRC, this.FRAG_SRC);
	if (!this.shaderLoc) {
    console.log(this.constructor.name + 
    						'.init() failed to create executable Shaders on the GPU. Bye!');
    return;
  }
// CUTE TRICK: let's print the NAME of this VBObox object: tells us which one!
//  else{console.log('You called: '+ this.constructor.name + '.init() fcn!');}

	gl.program = this.shaderLoc;		// (to match cuon-utils.js -- initShaders())

// b) Create VBO on GPU, fill it------------------------------------------------
	this.vboLoc = gl.createBuffer();	
  if (!this.vboLoc) {
    console.log(this.constructor.name + 
    						'.init() failed to create VBO in GPU. Bye!'); 
    return;
  }
  // Specify the purpose of our newly-created VBO on the GPU.  Your choices are:
  //	== "gl.ARRAY_BUFFER" : the VBO holds vertices, each made of attributes 
  // (positions, colors, normals, etc), or 
  //	== "gl.ELEMENT_ARRAY_BUFFER" : the VBO holds indices only; integer values 
  // that each select one vertex from a vertex array stored in another VBO.
  gl.bindBuffer(gl.ARRAY_BUFFER,	      // GLenum 'target' for this GPU buffer 
  								this.vboLoc);				  // the ID# the GPU uses for this buffer.

  // Fill the GPU's newly-created VBO object with the vertex data we stored in
  //  our 'vboContents' member (JavaScript Float32Array object).
  //  (Recall gl.bufferData() will evoke GPU's memory allocation & management: 
  //    use gl.bufferSubData() to modify VBO contents without changing VBO size)
  gl.bufferData(gl.ARRAY_BUFFER, 			  // GLenum target(same as 'bindBuffer()')
 					 				this.vboContents, 		// JavaScript Float32Array
  							 	gl.STATIC_DRAW);			// Usage hint.
  //	The 'hint' helps GPU allocate its shared memory for best speed & efficiency
  //	(see OpenGL ES specification for more info).  Your choices are:
  //		--STATIC_DRAW is for vertex buffers rendered many times, but whose 
  //				contents rarely or never change.
  //		--DYNAMIC_DRAW is for vertex buffers rendered many times, but whose 
  //				contents may change often as our program runs.
  //		--STREAM_DRAW is for vertex buffers that are rendered a small number of 
  // 			times and then discarded; for rapidly supplied & consumed VBOs.

  // c1) Find All Attributes:---------------------------------------------------
  //  Find & save the GPU location of all our shaders' attribute-variables and 
  //  uniform-variables (for switchToMe(), adjust(), draw(), reload(),etc.)
  this.a_PosLoc = gl.getAttribLocation(this.shaderLoc, 'a_Pos0');
  if(this.a_PosLoc < 0) {
    console.log(this.constructor.name + 
    						'.init() Failed to get GPU location of attribute a_Pos0');
    return -1;	// error exit.
  }
 	this.a_ColrLoc = gl.getAttribLocation(this.shaderLoc, 'a_Colr0');
  if(this.a_ColrLoc < 0) {
    console.log(this.constructor.name + 
    						'.init() failed to get the GPU location of attribute a_Colr0');
    return -1;	// error exit.
  }
  
  // c2) Find All Uniforms:-----------------------------------------------------
  //Get GPU storage location for each uniform var used in our shader programs: 
	this.u_ModelMatLoc = gl.getUniformLocation(this.shaderLoc, 'u_ModelMat0');
  if (!this.u_ModelMatLoc) { 
    console.log(this.constructor.name + 
    						'.init() failed to get GPU location for u_ModelMat1 uniform');
    return;
  }  
}

VBObox0.prototype.switchToMe = function() {
//==============================================================================
// Set GPU to use this VBObox's contents (VBO, shader, attributes, uniforms...)
//
// We only do this AFTER we called the init() function, which does the one-time-
// only setup tasks to put our VBObox contents into GPU memory.  !SURPRISE!
// even then, you are STILL not ready to draw our VBObox's contents onscreen!
// We must also first complete these steps:
//  a) tell the GPU to use our VBObox's shader program (already in GPU memory),
//  b) tell the GPU to use our VBObox's VBO  (already in GPU memory),
//  c) tell the GPU to connect the shader program's attributes to that VBO.

// a) select our shader program:
  gl.useProgram(this.shaderLoc);	
//		Each call to useProgram() selects a shader program from the GPU memory,
// but that's all -- it does nothing else!  Any previously used shader program's 
// connections to attributes and uniforms are now invalid, and thus we must now
// establish new connections between our shader program's attributes and the VBO
// we wish to use.  
  
// b) call bindBuffer to disconnect the GPU from its currently-bound VBO and
//  instead connect to our own already-created-&-filled VBO.  This new VBO can 
//    supply values to use as attributes in our newly-selected shader program:
	gl.bindBuffer(gl.ARRAY_BUFFER,	        // GLenum 'target' for this GPU buffer 
										this.vboLoc);			    // the ID# the GPU uses for our VBO.

// c) connect our newly-bound VBO to supply attribute variable values for each
// vertex to our SIMD shader program, using 'vertexAttribPointer()' function.
// this sets up data paths from VBO to our shader units:
  // 	Here's how to use the almost-identical OpenGL version of this function:
	//		http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribPointer.xml )
  gl.vertexAttribPointer(
		this.a_PosLoc,//index == ID# for the attribute var in your GLSL shader pgm;
		this.vboFcount_a_Pos0,// # of floats used by this attribute: 1,2,3 or 4?
		gl.FLOAT,			// type == what data type did we use for those numbers?
		false,				// isNormalized == are these fixed-point values that we need
									//									normalize before use? true or false
		this.vboStride,// Stride == #bytes we must skip in the VBO to move from the
		              // stored attrib for this vertex to the same stored attrib
		              //  for the next vertex in our VBO.  This is usually the 
									// number of bytes used to store one complete vertex.  If set 
									// to zero, the GPU gets attribute values sequentially from 
									// VBO, starting at 'Offset'.	
									// (Our vertex size in bytes: 4 floats for pos + 3 for color)
		this.vboOffset_a_Pos0);						
		              // Offset == how many bytes from START of buffer to the first
  								// value we will actually use?  (We start with position).
  gl.vertexAttribPointer(this.a_ColrLoc, this.vboFcount_a_Colr0, 
                        gl.FLOAT, false, 
                        this.vboStride, this.vboOffset_a_Colr0);
  							
// --Enable this assignment of each of these attributes to its' VBO source:
  gl.enableVertexAttribArray(this.a_PosLoc);
  gl.enableVertexAttribArray(this.a_ColrLoc);
}

VBObox0.prototype.isReady = function() {
//==============================================================================
// Returns 'true' if our WebGL rendering context ('gl') is ready to render using
// this objects VBO and shader program; else return false.
// see: https://developer.mozilla.org/en-US/docs/Web/API/WebGLRenderingContext/getParameter

var isOK = true;

  if(gl.getParameter(gl.CURRENT_PROGRAM) != this.shaderLoc)  {
    console.log(this.constructor.name + 
    						'.isReady() false: shader program at this.shaderLoc not in use!');
    isOK = false;
  }
  if(gl.getParameter(gl.ARRAY_BUFFER_BINDING) != this.vboLoc) {
      console.log(this.constructor.name + 
  						'.isReady() false: vbo at this.vboLoc not in use!');
    isOK = false;
  }
  return isOK;
}

VBObox0.prototype.adjust = function() {
//==============================================================================
// Update the GPU to newer, current values we now store for 'uniform' vars on 
// the GPU; and (if needed) update each attribute's stride and offset in VBO.

  // check: was WebGL context set to use our VBO & shader program?
  if(this.isReady()==false) {
        console.log('ERROR! before' + this.constructor.name + 
  						'.adjust() call you needed to call this.switchToMe()!!');
  }  
	// Adjust values for our uniforms,

		this.ModelMat.setIdentity();
// THIS DOESN'T WORK!!  this.ModelMatrix = g_worldMat;
  this.ModelMat.set(g_worldMat);	// use our global, shared camera.
// READY to draw in 'world' coord axes.
	
//  this.ModelMat.rotate(g_angleNow0, 0, 0, 1);	  // rotate drawing axes,
//  this.ModelMat.translate(0.35, 0, 0);							// then translate them.
  //  Transfer new uniforms' values to the GPU:-------------
  // Send  new 'ModelMat' values to the GPU's 'u_ModelMat1' uniform: 
  gl.uniformMatrix4fv(this.u_ModelMatLoc,	// GPU location of the uniform
  										false, 				// use matrix transpose instead?
  										this.ModelMat.elements);	// send data from Javascript.
  // Adjust the attributes' stride and offset (if necessary)
  // (use gl.vertexAttribPointer() calls and gl.enableVertexAttribArray() calls)
}

VBObox0.prototype.draw = function() {
//=============================================================================
// Render current VBObox contents.

  // check: was WebGL context set to use our VBO & shader program?
  if(this.isReady()==false) {
        console.log('ERROR! before' + this.constructor.name + 
  						'.draw() call you needed to call this.switchToMe()!!');
  }  
  // ----------------------------Draw the contents of the currently-bound VBO:
  gl.drawArrays(gl.LINES,				// use this drawing primitive, and
  0,	// start at this vertex number, and 
	this.vboVerts);	// draw this many vertices.
  }

VBObox0.prototype.reload = function() {
//=============================================================================
// Over-write current values in the GPU inside our already-created VBO: use 
// gl.bufferSubData() call to re-transfer some or all of our Float32Array 
// contents to our VBO without changing any GPU memory allocations.

 gl.bufferSubData(gl.ARRAY_BUFFER, 	// GLenum target(same as 'bindBuffer()')
                  0,                  // byte offset to where data replacement
                                      // begins in the VBO.
 					 				this.vboContents);   // the JS source-data array used to fill VBO

}
/*
VBObox0.prototype.empty = function() {
//=============================================================================
// Remove/release all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  However, make sure this step is reversible by a call to 
// 'restoreMe()': be sure to retain all our Float32Array data, all values for 
// uniforms, all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}

VBObox0.prototype.restore = function() {
//=============================================================================
// Replace/restore all GPU resources used by this VBObox object, including any 
// shader programs, attributes, uniforms, textures, samplers or other claims on 
// GPU memory.  Use our retained Float32Array data, all values for  uniforms, 
// all stride and offset values, etc.
//
//
// 		********   YOU WRITE THIS! ********
//
//
//
}
*/

const glsl = x => x;
//=============================================================================
//=============================================================================
function VBObox1() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox1' object that holds all data and fcns
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate shader program (a vertex-shader & fragment-shader pair) and one
// set of 'uniform' variables.

// Constructor goal: 
// Create and set member vars that will ELIMINATE ALL LITERALS (numerical values 
// written into code) in all other VBObox functions. Keeping all these (initial)
// values here, in this one coonstrutor function, ensures we can change them 
// easily WITHOUT disrupting any other code, ever!
  
	this.VERT_SRC =	//--------------------- VERTEX SHADER source code
  glsl`
  precision highp float;
  precision highp int;

  // materials
  uniform vec3 u_Ke;  //-- emission
  uniform vec3 u_Ka;  //-- ambient
  uniform vec3 u_Kd;  //-- diffuse
  uniform vec3 u_Ks;  //-- specular
  uniform int u_Kshiny; //-- exponent

  // light source
  uniform vec3 u_pos; //-- light position
  uniform vec3 u_Ia; //-- ambient light
  uniform vec3 u_Id; //-- diffuse light
  uniform vec3 u_Is; //-- specular light

  uniform mat4 u_ModelMatrix;
  uniform mat4 u_NormalMatrix;
  uniform mat4 u_MvpMatrix;
  uniform vec3 u_eyeWordPos;
  uniform bool u_isBlinn;
  uniform int u_fatt;

  attribute vec4 a_Pos1;
  attribute vec4 a_Norm1;

  varying vec4 v_Colr1;

  void main() {
    gl_Position = u_MvpMatrix * a_Pos1;
    vec4 v_Pos1 = u_ModelMatrix * a_Pos1;

    vec3 v_Norm1 = normalize(vec3(u_NormalMatrix * a_Norm1));
    vec3 N = normalize(v_Norm1);
    vec3 L = normalize(u_pos-vec3(v_Pos1));
    vec3 V = normalize(u_eyeWordPos - vec3(v_Pos1));
    float specular = 0.0;

    // attenuation
    float attenuation;
    float distance = length(u_pos-vec3(v_Pos1));

    if (u_fatt==1){
      attenuation = 1.0;
    }else if( u_fatt==2){
      attenuation = 1.0/(0.2*distance);
    }else if (u_fatt==3){
      attenuation = 1.0 / (1.0+0.1*distance+0.01*distance * distance);
    }

    // float diffuseCosine = max(dot(L, N), 0.0);
    float diffuseCosine = max(dot(L, N), 0.0)*attenuation;

    if (u_isBlinn){
      vec3 H = normalize(L + V); // Calculate halfway vector
      float specAngle=max(dot(H, N), 0.0);
      specular = pow(specAngle, float(u_Kshiny)); // Compute specular reflection
    }else{
      vec3 R = reflect(-L, N); // Calculate reflection vector
      float specAngle = max(dot(R, V), 0.0);
      specular = pow(specAngle, float(u_Kshiny)); // Compute specular reflection
  }
    vec3 ambient = u_Ia * u_Ka;
    vec3 emissive = u_Ke;
    vec3 diffuse = u_Id *diffuseCosine * u_Kd;
    vec3 specularHighlights = u_Is * specular * u_Ks;
    
    v_Colr1  = vec4(ambient + emissive + diffuse + specularHighlights, 1.0);
}

`;

//========YOUR CHOICE OF 3 Fragment shader programs=======
//				(use /* and */ to uncomment ONLY ONE)
// Each is an example of how to use the built-in vars for gl.POINTS to
// improve their on-screen appearance.
// a)'SQUARE points' -- DEFAULT; simple fixed-color square set by point-size.
// b) 'ROUND FLAT' -- uses 'gl_PointCoord' to make solid-color dot instead;
// c) 'SHADED Sphere' -- radial distance sets color to 'fake' a lit 3D sphere.
//   You too can be a 'shader writer'! What other fragment shaders would help?


 // a) SQUARE points:
	this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code 
  glsl`
  #ifdef GL_ES 
  precision mediump float;
  #endif
  varying vec4 v_Colr1;
  void main() {
    gl_FragColor = v_Colr1;
  }`;

  makeSphere();
  makeCylinder(0.1, 0.1);
  makeCube();
  makeCone();


  var mySize = sphVerts.length+cubeVerts.length+coneVerts.length

  this.vboContents = new Float32Array(mySize);
  
  sphStart = 0;
  for(i = 0, j = 0; j < sphVerts.length; i++, j++) {
    this.vboContents[i] = sphVerts[j];
  }
  coneStart = i;
  for(j=0; j< coneVerts.length; i++,j++) {
    this.vboContents[i] = coneVerts[j];
    }
  cubeStart = i;
  for (j = 0; j < cubeVerts.length; i++, j++) {
    this.vboContents[i] = cubeVerts[j];
  }

	this.vboVerts =this.vboContents.length / floatsPerVertex;	// # of vertices held in 'vboContents' array;
	this.FSIZE = this.vboContents.BYTES_PER_ELEMENT;  
	                              // bytes req'd by 1 vboContents array element;
																// (why? used to compute stride and offset 
																// in bytes for vertexAttribPointer() calls)
  this.vboBytes = this.vboContents.length * this.FSIZE;               
                                // (#  of floats in vboContents array) * 
                                // (# of bytes/float).
	this.vboStride = this.vboBytes / this.vboVerts;     
	                              // (== # of bytes to store one complete vertex).
	                              // From any attrib in a given vertex in the VBO, 
	                              // move forward by 'vboStride' bytes to arrive 
	                              // at the same attrib for the next vertex.
	                               
	            //----------------------Attribute sizes
  this.vboFcount_a_Pos1 =  3;    // # of floats in the VBO needed to store the
  // attribute named a_Pos1. (4: x,y,z values)
  this.vboFcount_a_Norm1 = 3;
  
  console.assert((this.vboFcount_a_Pos1 +     // check the size of each and
                  this.vboFcount_a_Norm1) *   // every attribute in our VBO
                  this.FSIZE == this.vboStride, // for agreeement with'stride'
                  "Uh oh! VBObox1.vboStride disagrees with attribute-size values!");
                  
  //             //----------------------Attribute offsets
	this.vboOffset_a_Pos1 = 0;    //# of bytes from START of vbo to the START
	                              // of 1st a_Pos1 attrib value in vboContents[]
  this.vboOffset_a_Norm1 =(this.vboOffset_a_Pos1) * this.FSIZE;    //# of bytes from START of vbo to the START
	                              // of 1st a_Pos1 attrib value in vboContents[]
  
  //-----------------------GPU memory locations:
	this.vboLoc;									// GPU Location for Vertex Buffer Object, 
	                              // returned by gl.createBuffer() function call
	this.shaderLoc;								// GPU Location for compiled Shader-program
	                            	// set by compile/link of VERT_SRC and FRAG_SRC.
	
  //------Attribute locations in our shaders:
	this.a_Pos1Loc;							  // GPU location: shader 'a_Pos1' attribute
	this.a_Norm1Loc;							// GPU location: shader 'a_Norm1' attribute

  this.ModelMatrix = new Matrix4();	// Transforms CVV axes to model axes.
  this.MvpMatrix = new Matrix4();	
  this.NormalMatrix = new Matrix4();	
  this.eyeWordPos = new Float32Array(3);

  this.u_ModelMatrixLoc;						// GPU location for u_ModelMat uniform
  this.u_NormalMatrixLoc;           
  this.u_MvpMatrixLoc; 
  this.u_eyeWordPosLoc;
  this.u_isBlinnLoc;
  this.u_fatt;

  this.setMaterial=function(matl) {
    gl.uniform3f(this.u_KeLoc,...matl.K_emit.slice(0,3));
    gl.uniform3f(this.u_KaLoc, ...matl.K_ambi.slice(0,3));
    gl.uniform3f(this.u_KdLoc, ...matl.K_diff.slice(0,3));
    gl.uniform3f(this.u_KsLoc, ...matl.K_spec.slice(0,3));
    gl.uniform1i(this.u_KshinyLoc, parseInt(g_shiny, 10));
    }

};


VBObox1.prototype.init = function() {
//==============================================================================
// Prepare the GPU to use all vertices, GLSL shaders, attributes, & uniforms
// kept in this VBObox. (This function usually called only once, within main()).
// Specifically:
// a) Create, compile, link our GLSL vertex- and fragment-shaders to form an
//  executable 'program' stored and ready to use inside the GPU.
// b) create a new VBO object in GPU memory and fill it by transferring in all
//  the vertex data held in our Float32array member 'VBOcontents'.
// c) Find & save the GPU location of all our shaders' attribute-variables and
//  uniform-variables (needed by switchToMe(), adjust(), draw(), reload(), etc.)
// -------------------
// CAREFUL!  before you can draw pictures using this VBObox contents,
//  you must call this VBObox object's switchToMe() function too!
//--------------------

  // a) Compile,link,upload shaders-----------------------------------------------
	this.shaderLoc = createProgram(gl, this.VERT_SRC, this.FRAG_SRC);
	if (!this.shaderLoc) {
    console.log(this.constructor.name + 
    						'.init() failed to create executable Shaders on the GPU. Bye!');
    return;
  }

// CUTE TRICK: let's print the NAME of this VBObox object: tells us which one!
//  else{console.log('You called: '+ this.constructor.name + '.init() fcn!');}

	gl.program = this.shaderLoc;		// (to match cuon-utils.js -- initShaders())

// b) Create VBO on GPU, fill it------------------------------------------------
	this.vboLoc = gl.createBuffer();	
  if (!this.vboLoc) {
    console.log(this.constructor.name + 
    						'.init() failed to create VBO in GPU. Bye!'); 
    return;
  }

  // c1) Find All Attributes:-----------------------------------------------------
  this.a_Pos1Loc = gl.getAttribLocation(this.shaderLoc, 'a_Pos1');
  this.a_Norm1Loc = gl.getAttribLocation(this.shaderLoc, 'a_Norm1');
  
  gl.bindBuffer(gl.ARRAY_BUFFER,	      // GLenum 'target' for this GPU buffer 
  								this.vboLoc);				  // the ID# the GPU uses for this buffer.
  											
  gl.bufferData(gl.ARRAY_BUFFER, 			  // GLenum target(same as 'bindBuffer()')
 					 				this.vboContents, 		// JavaScript Float32Array
    gl.STATIC_DRAW);			// Usage hint.
  

  // c2) Find All Uniforms:-----------------------------------------------------
  //Get GPU storage location for each uniform var used in our shader programs: 
  this.u_ModelMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_ModelMatrix');
  this.u_eyeWordPosLoc = gl.getUniformLocation(this.shaderLoc, 'u_eyeWordPos');
  this.u_MvpMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_MvpMatrix');
  this.u_NormalMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_NormalMatrix');
  this.u_isBlinnLoc = gl.getUniformLocation(this.shaderLoc, 'u_isBlinn');
  this.u_fatt = gl.getUniformLocation(this.shaderLoc, 'u_fatt');

  // console.log(this.u_fattLoc)
  
  if (!this.u_ModelMatrixLoc || !this.u_eyeWordPosLoc || !this.u_MvpMatrixLoc
     ||!this.u_NormalMatrixLoc || !this.u_isBlinnLoc || !this.u_fatt) { 
    console.log(this.constructor.name + 
    						'.init() failed to get GPU location for uniform locations');
    return;
  }


  this.u_KeLoc = gl.getUniformLocation(this.shaderLoc, 'u_Ke');
  this.u_KaLoc = gl.getUniformLocation(this.shaderLoc, 'u_Ka');
  this.u_KdLoc = gl.getUniformLocation(this.shaderLoc, 'u_Kd');
  this.u_KsLoc = gl.getUniformLocation(this.shaderLoc, 'u_Ks');
  this.u_KshinyLoc = gl.getUniformLocation(this.shaderLoc, 'u_Kshiny');
  
  if (!this.u_KeLoc || !this.u_KaLoc || !this.u_KdLoc || !this.u_KshinyLoc) {
    console.log('Failed to get one or more material storage locations');
    return;
  }  

  this.u_pos = gl.getUniformLocation(this.shaderLoc,  'u_pos');
  this.u_diff = gl.getUniformLocation(this.shaderLoc, 'u_Id');
  this.u_ambi = gl.getUniformLocation(this.shaderLoc, 'u_Ia');
  this.u_spec = gl.getUniformLocation(this.shaderLoc, 'u_Is');

  if( !this.u_pos || !this.u_diff || !this.u_ambi || !this.u_spec) {
      console.log(this.constructor.name + ' failed to get one or more lighting uniform storage locations.');
      return;
  }


}


VBObox1.prototype.switchToMe = function () {
//==============================================================================
// Set GPU to use this VBObox's contents (VBO, shader, attributes, uniforms...)
//
// We only do this AFTER we called the init() function, which does the one-time-
// only setup tasks to put our VBObox contents into GPU memory.  !SURPRISE!
// even then, you are STILL not ready to draw our VBObox's contents onscreen!
// We must also first complete these steps:
//  a) tell the GPU to use our VBObox's shader program (already in GPU memory),
//  b) tell the GPU to use our VBObox's VBO  (already in GPU memory),
//  c) tell the GPU to connect the shader program's attributes to that VBO.

// a) select our shader program:
gl.useProgram(this.shaderLoc);	
//		Each call to useProgram() selects a shader program from the GPU memory,
// but that's all -- it does nothing else!  Any previously used shader program's
// connections to attributes and uniforms are now invalid, and thus we must now
// establish new connections between our shader program's attributes and the VBO
// we wish to use.

gl.uniform1i(this.u_isBlinnLoc, g_isBlinn);
gl.uniform1i(this.u_fatt, g_fatt);
gl.uniform3fv(this.u_eyeWordPosLoc, this.eyeWordPos);

  
this.setMaterial(matl1)

// light source
gl.uniform3f(this.u_pos, g_lightPosX, g_lightPosY, g_lightPosZ);
gl.uniform3f(this.u_diff, g_lightDiffR, g_lightDiffG, g_lightDiffB);
gl.uniform3f(this.u_ambi, g_lightAmbiR, g_lightAmbiG, g_lightAmbiB);
gl.uniform3f(this.u_spec, g_lightSpecR, g_lightSpecG, g_lightSpecB);
  
// b) call bindBuffer to disconnect the GPU from its currently-bound VBO and
//  instead connect to our own already-created-&-filled VBO.  This new VBO can 
//    supply values to use as attributes in our newly-selected shader program:
	gl.bindBuffer(gl.ARRAY_BUFFER,	    // GLenum 'target' for this GPU buffer 
										this.vboLoc);			// the ID# the GPU uses for our VBO.

  
  gl.vertexAttribPointer(
    this.a_Pos1Loc,
    this.vboFcount_a_Pos1, 
    gl.FLOAT,
    false,
    this.vboStride,
    this.vboOffset_a_Pos1);	
  
  gl.vertexAttribPointer(
    this.a_Norm1Loc,
    this.vboFcount_a_Norm1,             
    gl.FLOAT,
    false, 
    this.vboStride,
    this.vboOffset_a_Norm1);	
  
  gl.enableVertexAttribArray(this.a_Pos1Loc);
  gl.enableVertexAttribArray(this.a_Norm1Loc);
}

VBObox1.prototype.isReady = function() {
//==============================================================================
// Returns 'true' if our WebGL rendering context ('gl') is ready to render using
// this objects VBO and shader program; else return false.
// see: https://developer.mozilla.org/en-US/docs/Web/API/WebGLRenderingContext/getParameter

var isOK = true;

  if(gl.getParameter(gl.CURRENT_PROGRAM) != this.shaderLoc)  {
    console.log(this.constructor.name + 
    						'.isReady() false: shader program at this.shaderLoc not in use!');
    isOK = false;
  }
  if(gl.getParameter(gl.ARRAY_BUFFER_BINDING) != this.vboLoc) {
      console.log(this.constructor.name + 
  						'.isReady() false: vbo at this.vboLoc not in use!');
    isOK = false;
  }
  return isOK;
}

VBObox1.prototype.adjust = function() {
  //==============================================================================
  // Update the GPU to newer, current values we now store for 'uniform' vars on 
  // the GPU; and (if needed) update each attribute's stride and offset in VBO.
  
    // check: was WebGL context set to use our VBO & shader program?
    if(this.isReady()==false) {
          console.log('ERROR! before' + this.constructor.name + 
                '.adjust() call you needed to call this.switchToMe()!!');
    }
    // Adjust values for our uniforms,
    this.eyeWordPos.set([g_EyeX, g_EyeY, g_EyeZ]);
    this.ModelMatrix.setIdentity();
    this.MvpMatrix.setIdentity();
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.rotate(g_angleNow1, 0, 0, 1);	// spin drawing axes,
  
  //  this.ModelMatrix.rotate(g_angleNow1, 0, 0, 1);	// -spin drawing axes,
  // this.ModelMatrix.translate(1.0, -2.0, 0); // then translate them.
    
  drawSegment(
    gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
    this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  
  // shape 2
  this.MvpMatrix.set(g_worldMat);
  this.ModelMatrix.setIdentity();
  
  this.ModelMatrix.translate(-2.0, 0.8, 0.0);
  this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
  this.ModelMatrix.rotate(g_angleNow1, 0, 1, 0);	// spin drawing axes,
  pushMatrix(this.ModelMatrix);
    //ModelMatrix.rotate(currentAngle, 1, 0, 0);  // spin around y axis.
    this.ModelMatrix.scale(0.2, 0.7, 0.7);
    this.setMaterial(matl1)
    drawSegment(
      gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  this.ModelMatrix = popMatrix();
    
  pushMatrix(this.ModelMatrix);
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.translate(0.0, 2.0, 0.0)
    this.ModelMatrix.rotate(-90.0, 1.0, 0.0, 0.0)
    drawSegment(
      gl, coneStart, coneVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  this.ModelMatrix = popMatrix(); 
    
  
  // flex shape  1
  this.MvpMatrix.set(g_worldMat);
  this.ModelMatrix.setIdentity();
  
  this.ModelMatrix.translate(-2.0, 3.5, 0.0);
  this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
  this.ModelMatrix.rotate(g_angleNow1, 0, 1, 0);	// spin drawing axes,
  pushMatrix(this.ModelMatrix);
    //ModelMatrix.rotate(currentAngle, 1, 0, 0);  // spin around y axis.
    this.ModelMatrix.scale(0.2, 0.7, 0.7);
    this.setMaterial(matl1)
    drawSegment(
      gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  this.ModelMatrix = popMatrix();
    
  pushMatrix(this.ModelMatrix);
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.translate(0.0, 2.0, 0.0)
  this.ModelMatrix.rotate(-90.0, 1.0, 0.0, 0.0)
    drawSegment(
      gl, coneStart, coneVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
    
    // top
    this.setMaterial(matl1)
    this.MvpMatrix.set(g_worldMat);

  this.ModelMatrix.scale(1.5, 1.5, 0.5);
  this.ModelMatrix.translate(0.0, -1.0, -2.0)
  this.ModelMatrix.rotate(angleA, 0.0, 1.0, 0.0)
    drawSegment(
        gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
        this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  
  this.ModelMatrix = popMatrix(); 
    

  // flex shape 2
  this.MvpMatrix.set(g_worldMat);
  this.ModelMatrix.setIdentity();
  
  this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
  this.ModelMatrix.translate(1.0, 0.8, -3.0);
    this.ModelMatrix.scale(0.5, 0.5, 0.5);
  
  this.setMaterial(matl2)
    this.ModelMatrix.rotate(angleA, 1, 0, 0);  // spin around y axis.
    this.ModelMatrix.scale(0.2, 0.7, 0.7);
    drawSegment(
      gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
    
  pushMatrix(this.ModelMatrix);
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.translate(0.0, 2.0, 0.0)
    this.ModelMatrix.rotate(angleB, 1.0, 0.0, 0.0)
    drawSegment(
      gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)

  this.MvpMatrix.set(g_worldMat);
  this.ModelMatrix.translate(0.0, 2.0, 0.0)
  this.ModelMatrix.rotate(angleC, 1.0, 0.0, 0.0)
  drawSegment(
    gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
    this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  
  // top
  this.setMaterial(matl1)
  this.MvpMatrix.set(g_worldMat);
  this.ModelMatrix.translate(0.0, 2.0, 0.0)
  drawSegment(
    gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
    this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)

  this.ModelMatrix = popMatrix(); 
    

  // flex shape 3

    // flex shape 2
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.setIdentity();
    
    this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
    this.ModelMatrix.translate(3.0, 0.8, -3.0);
    this.ModelMatrix.scale(1.5, 0.5, 0.5);
    
    this.setMaterial(matl2)
      this.ModelMatrix.rotate(angleA, 1, 0, 0);  // spin around y axis.
      this.ModelMatrix.scale(0.2, 0.7, 0.7);
      drawSegment(
        gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
        this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
    
    pushMatrix(this.ModelMatrix);
      this.MvpMatrix.set(g_worldMat);
      this.ModelMatrix.translate(0.0, 2.0, 0.0)
      this.ModelMatrix.rotate(angleB, 0.0, 0.0, 1.0)
      drawSegment(
        gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
        this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
    
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.translate(0.0, 2.0, 0.0)
    this.ModelMatrix.rotate(angleC, 1.0, 1.0, 0.0)
    drawSegment(
      gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  
    // top
    this.MvpMatrix.set(g_worldMat);
    this.ModelMatrix.translate(0.0, 2.0, 0.0)
    drawSegment(
      gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
      this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
  
    this.ModelMatrix = popMatrix(); 
      
    
  }

function drawSegment(gl, start, Vertices, ModelMatrix, MvpMatrix, NormalMatrix, u_NormalMatrixLoc, u_MvpMatrixLoc, u_ModelMatrixLoc) {
  NormalMatrix.setInverseOf(ModelMatrix);
  NormalMatrix.transpose();

  MvpMatrix.multiply(ModelMatrix);
  gl.uniformMatrix4fv(u_ModelMatrixLoc,	false, ModelMatrix.elements);									
  gl.uniformMatrix4fv(u_NormalMatrixLoc, false,NormalMatrix.elements);
  gl.uniformMatrix4fv(u_MvpMatrixLoc, false, MvpMatrix.elements);
    
  // Draw just the the cylinder's vertices:
  return gl.drawArrays(gl.TRIANGLE_STRIP,				// use this drawing primitive, and
                start/floatsPerVertex, // start at this vertex number, and
    Vertices.length / floatsPerVertex);	// draw this many vertices.

}



VBObox1.prototype.draw = function() {
//=============================================================================
// Send commands to GPU to select and render current VBObox contents.
  // check: was WebGL context set to use our VBO & shader program?
  if(this.isReady()==false) {
        console.log('ERROR! before' + this.constructor.name + 
  						'.draw() call you needed to call this.switchToMe()!!');
  }

  // ----------------------------Draw the contents of the currently-bound VBO:
  //gl.drawArrays(gl.TRIANGLE_STRIP, sphStart / floatsPerVertex, sphVerts.length / floatsPerVertex);
}


VBObox1.prototype.reload = function() {
//=============================================================================
// Over-write current values in the GPU for our already-created VBO: use 
// gl.bufferSubData() call to re-transfer some or all of our Float32Array 
// contents to our VBO without changing any GPU memory allocations.

 gl.bufferSubData(gl.ARRAY_BUFFER, 	// GLenum target(same as 'bindBuffer()')
                  0,                  // byte offset to where data replacement
                                      // begins in the VBO.
 					 				this.vboContents);   // the JS source-data array used to fill VBO
}


//=============================================================================
//=============================================================================
function VBObox2() {
//=============================================================================
//=============================================================================
// CONSTRUCTOR for one re-usable 'VBObox2' object that holds all data and fcns
// needed to render vertices from one Vertex Buffer Object (VBO) using one 
// separate shader program (a vertex-shader & fragment-shader pair) and one
// set of 'uniform' variables.

// Constructor goal: 
// Create and set member vars that will ELIMINATE ALL LITERALS (numerical values 
// written into code) in all other VBObox functions. Keeping all these (initial)
// values here, in this one coonstrutor function, ensures we can change them 
// easily WITHOUT disrupting any other code, ever!
  
this.VERT_SRC =	//--------------------- VERTEX SHADER source code
glsl`
precision highp float;
precision highp int;


uniform vec3 u_Kd;  //-- diffuse
uniform mat4 u_ModelMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_MvpMatrix;

attribute vec4 a_Pos2;
attribute vec4 a_Norm2;

varying vec3 v_Kd;
varying vec4 v_Pos2;
varying vec3 v_Norm2;

void main() {
  v_Kd=u_Kd;
  gl_Position = u_MvpMatrix * a_Pos2;
  v_Pos2 = u_ModelMatrix * a_Pos2;
  v_Norm2 = normalize(vec3(u_NormalMatrix * a_Norm2));
}
`;

	this.FRAG_SRC = //---------------------- FRAGMENT SHADER source code 
  glsl`
  #ifdef GL_ES 
  precision highp float;
  precision highp int;
  #endif

  // materials
  uniform vec3 u_Ke;  //-- emission
  uniform vec3 u_Ka;  //-- ambient
  uniform vec3 u_Ks;  //-- specular
  uniform int u_Kshiny; //-- exponent

  // light source
  uniform vec3 u_pos;
  uniform vec3 u_Ia; //-- ambient light
  uniform vec3 u_Id; //-- diffuse light
  uniform vec3 u_Is; //-- specular light

  uniform vec3 u_eyeWordPos;
  uniform bool u_isBlinn;
  uniform int u_fatt;

  varying vec3 v_Kd;
  varying vec4 v_Pos2;
  varying vec3 v_Norm2;

  void main() {
    vec3 N = normalize(v_Norm2);
    vec3 L = normalize(u_pos-vec3(v_Pos2));
    vec3 V = normalize(u_eyeWordPos - vec3(v_Pos2));
    float specular = 0.0;

    // attenuation
    float attenuation;
    float distance = length(u_pos-vec3(v_Pos2));

    if (u_fatt==1){
      attenuation = 1.0;
    }else if( u_fatt==2){
      attenuation = 1.0/(0.2*distance);
    }else if (u_fatt==3){
      attenuation = 1.0 / (1.0+0.1*distance+0.01*distance * distance);
    }

    // float diffuseCosine = max(dot(L, N), 0.0);
    float diffuseCosine = max(dot(L, N), 0.0)*attenuation;

    if (u_isBlinn){
      vec3 H = normalize(L + V); // Calculate halfway vector
      float specAngle=max(dot(H, N), 0.0);
      specular = pow(specAngle, float(u_Kshiny)); // Compute specular reflection
    }else{
      vec3 R = reflect(-L, N); // Calculate reflection vector
      float specAngle = max(dot(R, V), 0.0);
      specular = pow(specAngle, float(u_Kshiny)); // Compute specular reflection
  }
    vec3 ambient = u_Ia * u_Ka;
    vec3 emissive = u_Ke;
    vec3 diffuse = u_Id *diffuseCosine * v_Kd;
    vec3 specularHighlights = u_Is * specular * u_Ks;

    gl_FragColor  = vec4(ambient + emissive + diffuse + specularHighlights, 1.0);
}
`;

  makeSphere();
  makeCylinder(0.1, 0.1);
  makeCube();
  makeCone();


  var mySize = sphVerts.length+cubeVerts.length+coneVerts.length

  this.vboContents = new Float32Array(mySize);

  sphStart = 0;
  for(i = 0, j = 0; j < sphVerts.length; i++, j++) {
    this.vboContents[i] = sphVerts[j];
  }
  coneStart = i;
  for(j=0; j< coneVerts.length; i++,j++) {
    this.vboContents[i] = coneVerts[j];
    }
  cubeStart = i;
  for (j = 0; j < cubeVerts.length; i++, j++) {
    this.vboContents[i] = cubeVerts[j];
  }

  this.vboVerts =this.vboContents.length / floatsPerVertex;	// # of vertices held in 'vboContents' array;

	this.FSIZE = this.vboContents.BYTES_PER_ELEMENT;
	                              // bytes req'd by 1 vboContents array element;
																// (why? used to compute stride and offset 
																// in bytes for vertexAttribPointer() calls)
  this.vboBytes = this.vboContents.length * this.FSIZE;               
                                // (#  of floats in vboContents array) * 
                                // (# of bytes/float).
	this.vboStride = this.vboBytes / this.vboVerts;     
	                              // From any attrib in a given vertex, 
	                              // move forward by 'vboStride' bytes to arrive 
	                              // at the same attrib for the next vertex. 
	                              // (== # of bytes used to store one vertex) 
	                              
	            //----------------------Attribute sizes
  this.vboFcount_a_Pos2 =  3;    // # of floats in the VBO needed to store the
                // attribute named a_Pos1. (4: x,y,z,w values)
  this.vboFcount_a_Norm2 = 3;  // # of floats for this attrib (just one!)   
  console.assert(((this.vboFcount_a_Pos2 +     // check the size of each and 
  this.vboFcount_a_Norm2) *   // every attribute in our VBO
  this.FSIZE == this.vboStride), // for agreeement with'stride'
  "Uh oh! VBObox1.vboStride disagrees with attribute-size values!");

	this.vboOffset_a_Pos2 = 0;   
	                              //# of bytes from START of vbo to the START
	                              // of 1st a_Position attrib value in vboContents[]
  this.vboOffset_a_Norm2 = (this.vboFcount_a_Pos2) * this.FSIZE;  
                                // == 4 floats * bytes/float
                                //# of bytes from START of vbo to the START
                                // of 1st a_Color attrib value in vboContents[]
  //-----------------------GPU memory locations:
	this.vboLoc;									// GPU Location for Vertex Buffer Object, 
	                              // returned by gl.createBuffer() function call
	this.shaderLoc;								// GPU Location for compiled Shader-program
	                            	// set by compile/link of VERT_SRC and FRAG_SRC.
	
  //------Attribute locations in our shaders:
	this.a_Pos2Loc;							  // GPU location: shader 'a_Pos1' attribute
	this.a_Norm2Loc;							// GPU location: shader 'a_Norm1' attribute

  this.ModelMatrix = new Matrix4();	// Transforms CVV axes to model axes.
  this.MvpMatrix = new Matrix4();	
  this.NormalMatrix = new Matrix4();	
  this.eyeWordPos = new Float32Array(3);

  this.u_ModelMatrixLoc;						// GPU location for u_ModelMat uniform
  this.u_NormalMatrixLoc;           
  this.u_MvpMatrixLoc; 
  this.u_eyeWordPosLoc;
  this.u_isBlinnLoc;
  this.u_fatt;

  this.setMaterial=function(matl) {
    gl.uniform3f(this.u_KeLoc,...matl.K_emit.slice(0,3));
    gl.uniform3f(this.u_KaLoc, ...matl.K_ambi.slice(0,3));
    gl.uniform3f(this.u_KdLoc, ...matl.K_diff.slice(0,3));
    gl.uniform3f(this.u_KsLoc, ...matl.K_spec.slice(0,3));
    gl.uniform1i(this.u_KshinyLoc, parseInt(g_shiny, 10));
    }
};


VBObox2.prototype.init = function() {
  //==============================================================================
  // Prepare the GPU to use all vertices, GLSL shaders, attributes, & uniforms
  // kept in this VBObox. (This function usually called only once, within main()).
  // Specifically:
  // a) Create, compile, link our GLSL vertex- and fragment-shaders to form an
  //  executable 'program' stored and ready to use inside the GPU.
  // b) create a new VBO object in GPU memory and fill it by transferring in all
  //  the vertex data held in our Float32array member 'VBOcontents'.
  // c) Find & save the GPU location of all our shaders' attribute-variables and
  //  uniform-variables (needed by switchToMe(), adjust(), draw(), reload(), etc.)
  // -------------------
  // CAREFUL!  before you can draw pictures using this VBObox contents,
  //  you must call this VBObox object's switchToMe() function too!
  //--------------------
  
    // a) Compile,link,upload shaders-----------------------------------------------
    this.shaderLoc = createProgram(gl, this.VERT_SRC, this.FRAG_SRC);
    if (!this.shaderLoc) {
      console.log(this.constructor.name + 
                  '.init() failed to create executable Shaders on the GPU. Bye!');
      return;
    }
  
  // CUTE TRICK: let's print the NAME of this VBObox object: tells us which one!
  //  else{console.log('You called: '+ this.constructor.name + '.init() fcn!');}
  
    gl.program = this.shaderLoc;		// (to match cuon-utils.js -- initShaders())
  
  // b) Create VBO on GPU, fill it------------------------------------------------
    this.vboLoc = gl.createBuffer();	
    if (!this.vboLoc) {
      console.log(this.constructor.name + 
                  '.init() failed to create VBO in GPU. Bye!'); 
      return;
    }
  
    // c1) Find All Attributes:-----------------------------------------------------
    this.a_Pos2Loc = gl.getAttribLocation(this.shaderLoc, 'a_Pos2');
    this.a_Norm2Loc = gl.getAttribLocation(this.shaderLoc, 'a_Norm2');
    
    gl.bindBuffer(gl.ARRAY_BUFFER,	      // GLenum 'target' for this GPU buffer 
                    this.vboLoc);				  // the ID# the GPU uses for this buffer.
                          
    gl.bufferData(gl.ARRAY_BUFFER, 			  // GLenum target(same as 'bindBuffer()')
                      this.vboContents, 		// JavaScript Float32Array
      gl.STATIC_DRAW);			// Usage hint.
    
  
    // c2) Find All Uniforms:-----------------------------------------------------
    //Get GPU storage location for each uniform var used in our shader programs: 
    this.u_ModelMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_ModelMatrix');
    this.u_eyeWordPosLoc = gl.getUniformLocation(this.shaderLoc, 'u_eyeWordPos');
    this.u_MvpMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_MvpMatrix');
    this.u_NormalMatrixLoc = gl.getUniformLocation(this.shaderLoc, 'u_NormalMatrix');
    this.u_isBlinnLoc = gl.getUniformLocation(this.shaderLoc, 'u_isBlinn');
    this.u_fatt = gl.getUniformLocation(this.shaderLoc, 'u_fatt');
  
    // console.log(this.u_fattLoc)
    
    if (!this.u_ModelMatrixLoc || !this.u_eyeWordPosLoc || !this.u_MvpMatrixLoc
       ||!this.u_NormalMatrixLoc || !this.u_isBlinnLoc || !this.u_fatt) { 
      console.log(this.constructor.name + 
                  '.init() failed to get GPU location for uniform locations');
      return;
    }
  
  
    this.u_KeLoc = gl.getUniformLocation(this.shaderLoc, 'u_Ke');
    this.u_KaLoc = gl.getUniformLocation(this.shaderLoc, 'u_Ka');
    this.u_KdLoc = gl.getUniformLocation(this.shaderLoc, 'u_Kd');
    this.u_KsLoc = gl.getUniformLocation(this.shaderLoc, 'u_Ks');
    this.u_KshinyLoc = gl.getUniformLocation(this.shaderLoc, 'u_Kshiny');
    
    if (!this.u_KeLoc || !this.u_KaLoc || !this.u_KdLoc || !this.u_KshinyLoc) {
      console.log('Failed to get one or more material storage locations');
      return;
    }  
  
    this.u_pos = gl.getUniformLocation(this.shaderLoc,  'u_pos');
    this.u_diff = gl.getUniformLocation(this.shaderLoc, 'u_Id');
    this.u_ambi = gl.getUniformLocation(this.shaderLoc, 'u_Ia');
    this.u_spec = gl.getUniformLocation(this.shaderLoc, 'u_Is');
  
    if( !this.u_pos || !this.u_diff || !this.u_ambi || !this.u_spec) {
        console.log(this.constructor.name + ' failed to get one or more lighting uniform storage locations.');
        return;
    }
  
  
  }

VBObox2.prototype.switchToMe = function () {
    //==============================================================================
    // Set GPU to use this VBObox's contents (VBO, shader, attributes, uniforms...)
    //
    // We only do this AFTER we called the init() function, which does the one-time-
    // only setup tasks to put our VBObox contents into GPU memory.  !SURPRISE!
    // even then, you are STILL not ready to draw our VBObox's contents onscreen!
    // We must also first complete these steps:
    //  a) tell the GPU to use our VBObox's shader program (already in GPU memory),
    //  b) tell the GPU to use our VBObox's VBO  (already in GPU memory),
    //  c) tell the GPU to connect the shader program's attributes to that VBO.
    
    // a) select our shader program:
    gl.useProgram(this.shaderLoc);	
    //		Each call to useProgram() selects a shader program from the GPU memory,
    // but that's all -- it does nothing else!  Any previously used shader program's
    // connections to attributes and uniforms are now invalid, and thus we must now
    // establish new connections between our shader program's attributes and the VBO
    // we wish to use.
      
    gl.uniform1i(this.u_isBlinnLoc, g_isBlinn);
    gl.uniform1i(this.u_fatt, g_fatt);
    gl.uniform3fv(this.u_eyeWordPosLoc, this.eyeWordPos);
    
      
    this.setMaterial(matl1)
    
    // light source
    gl.uniform3f(this.u_pos, g_lightPosX, g_lightPosY, g_lightPosZ);
    gl.uniform3f(this.u_diff, g_lightDiffR, g_lightDiffG, g_lightDiffB);
    gl.uniform3f(this.u_ambi, g_lightAmbiR, g_lightAmbiG, g_lightAmbiB);
    gl.uniform3f(this.u_spec, g_lightSpecR, g_lightSpecG, g_lightSpecB);
      
    // b) call bindBuffer to disconnect the GPU from its currently-bound VBO and
    //  instead connect to our own already-created-&-filled VBO.  This new VBO can 
    //    supply values to use as attributes in our newly-selected shader program:
      gl.bindBuffer(gl.ARRAY_BUFFER,	    // GLenum 'target' for this GPU buffer 
                        this.vboLoc);			// the ID# the GPU uses for our VBO.
    
      
      gl.vertexAttribPointer(
        this.a_Pos2Loc,
        this.vboFcount_a_Pos2, 
        gl.FLOAT,
        false,
        this.vboStride,
        this.vboOffset_a_Pos2);	
      
      gl.vertexAttribPointer(
        this.a_Norm2Loc,
        this.vboFcount_a_Norm2,             
        gl.FLOAT,
        false, 
        this.vboStride,
        this.vboOffset_a_Norm2);	
      
      gl.enableVertexAttribArray(this.a_Pos2Loc);
      gl.enableVertexAttribArray(this.a_Norm2Loc);
    }
    
VBObox2.prototype.isReady = function() {
    //==============================================================================
    // Returns 'true' if our WebGL rendering context ('gl') is ready to render using
    // this objects VBO and shader program; else return false.
    // see: https://developer.mozilla.org/en-US/docs/Web/API/WebGLRenderingContext/getParameter
    
    var isOK = true;
    
      if(gl.getParameter(gl.CURRENT_PROGRAM) != this.shaderLoc)  {
        console.log(this.constructor.name + 
                    '.isReady() false: shader program at this.shaderLoc not in use!');
        isOK = false;
      }
      if(gl.getParameter(gl.ARRAY_BUFFER_BINDING) != this.vboLoc) {
          console.log(this.constructor.name + 
                  '.isReady() false: vbo at this.vboLoc not in use!');
        isOK = false;
      }
      return isOK;
    }
    
VBObox2.prototype.adjust = function() {
      //==============================================================================
      // Update the GPU to newer, current values we now store for 'uniform' vars on 
      // the GPU; and (if needed) update each attribute's stride and offset in VBO.
      
        // check: was WebGL context set to use our VBO & shader program?
        if(this.isReady()==false) {
              console.log('ERROR! before' + this.constructor.name + 
                    '.adjust() call you needed to call this.switchToMe()!!');
        }
        // Adjust values for our uniforms,
        this.eyeWordPos.set([g_EyeX, g_EyeY, g_EyeZ]);
        this.ModelMatrix.setIdentity();
        this.MvpMatrix.setIdentity();
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.rotate(g_angleNow1, 0, 0, 1);	// spin drawing axes,
      
      //  this.ModelMatrix.rotate(g_angleNow1, 0, 0, 1);	// -spin drawing axes,
      // this.ModelMatrix.translate(1.0, -2.0, 0); // then translate them.
        
      drawSegment(
        gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
        this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      
      // shape 2
      this.MvpMatrix.set(g_worldMat);
      this.ModelMatrix.setIdentity();
      
      this.ModelMatrix.translate(-2.0, 0.8, 0.0);
      this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
      this.ModelMatrix.rotate(g_angleNow1, 0, 1, 0);	// spin drawing axes,
      pushMatrix(this.ModelMatrix);
        //ModelMatrix.rotate(currentAngle, 1, 0, 0);  // spin around y axis.
        this.ModelMatrix.scale(0.2, 0.7, 0.7);
        this.setMaterial(matl1)
        drawSegment(
          gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      this.ModelMatrix = popMatrix();
        
      pushMatrix(this.ModelMatrix);
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.translate(0.0, 2.0, 0.0)
        this.ModelMatrix.rotate(-90.0, 1.0, 0.0, 0.0)
        drawSegment(
          gl, coneStart, coneVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      this.ModelMatrix = popMatrix(); 
        
      
      // flex shape  1
      this.MvpMatrix.set(g_worldMat);
      this.ModelMatrix.setIdentity();
      
      this.ModelMatrix.translate(-2.0, 3.5, 0.0);
      this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
      this.ModelMatrix.rotate(g_angleNow1, 0, 1, 0);	// spin drawing axes,
      pushMatrix(this.ModelMatrix);
        //ModelMatrix.rotate(currentAngle, 1, 0, 0);  // spin around y axis.
        this.ModelMatrix.scale(0.2, 0.7, 0.7);
        this.setMaterial(matl1)
        drawSegment(
          gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      this.ModelMatrix = popMatrix();
        
      pushMatrix(this.ModelMatrix);
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.translate(0.0, 2.0, 0.0)
      this.ModelMatrix.rotate(-90.0, 1.0, 0.0, 0.0)
        drawSegment(
          gl, coneStart, coneVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
        
        // top
        this.setMaterial(matl1)
        this.MvpMatrix.set(g_worldMat);
    
      this.ModelMatrix.scale(1.5, 1.5, 0.5);
      this.ModelMatrix.translate(0.0, -1.0, -2.0)
      this.ModelMatrix.rotate(angleA, 0.0, 1.0, 0.0)
        drawSegment(
            gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
            this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      
      this.ModelMatrix = popMatrix(); 
        
    
      // flex shape 2
      this.MvpMatrix.set(g_worldMat);
      this.ModelMatrix.setIdentity();
      
      this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
      this.ModelMatrix.translate(1.0, 0.8, -3.0);
        this.ModelMatrix.scale(0.5, 0.5, 0.5);
      
      this.setMaterial(matl2)
        this.ModelMatrix.rotate(angleA, 1, 0, 0);  // spin around y axis.
        this.ModelMatrix.scale(0.2, 0.7, 0.7);
        drawSegment(
          gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
        
      pushMatrix(this.ModelMatrix);
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.translate(0.0, 2.0, 0.0)
        this.ModelMatrix.rotate(angleB, 1.0, 0.0, 0.0)
        drawSegment(
          gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
    
      this.MvpMatrix.set(g_worldMat);
      this.ModelMatrix.translate(0.0, 2.0, 0.0)
      this.ModelMatrix.rotate(angleC, 1.0, 0.0, 0.0)
      drawSegment(
        gl, cubeStart, cubeVerts, this.ModelMatrix, this.MvpMatrix,
        this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      
      // top
      this.setMaterial(matl1)
      this.MvpMatrix.set(g_worldMat);
      this.ModelMatrix.translate(0.0, 2.0, 0.0)
      drawSegment(
        gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
        this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
    
      this.ModelMatrix = popMatrix(); 
        
    
      // flex shape 3
    
        // flex shape 2
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.setIdentity();
        
        this.ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
        this.ModelMatrix.translate(3.0, 0.8, -3.0);
        this.ModelMatrix.scale(1.5, 0.5, 0.5);
        
        this.setMaterial(matl2)
          this.ModelMatrix.rotate(angleA, 1, 0, 0);  // spin around y axis.
          this.ModelMatrix.scale(0.2, 0.7, 0.7);
          drawSegment(
            gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
            this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
        
        pushMatrix(this.ModelMatrix);
          this.MvpMatrix.set(g_worldMat);
          this.ModelMatrix.translate(0.0, 2.0, 0.0)
          this.ModelMatrix.rotate(angleB, 0.0, 0.0, 1.0)
          drawSegment(
            gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
            this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
        
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.translate(0.0, 2.0, 0.0)
        this.ModelMatrix.rotate(angleC, 1.0, 1.0, 0.0)
        drawSegment(
          gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      
        // top
        this.MvpMatrix.set(g_worldMat);
        this.ModelMatrix.translate(0.0, 2.0, 0.0)
        drawSegment(
          gl, sphStart, sphVerts, this.ModelMatrix, this.MvpMatrix,
          this.NormalMatrix, this.u_NormalMatrixLoc, this.u_MvpMatrixLoc, this.u_ModelMatrixLoc)
      
        this.ModelMatrix = popMatrix(); 
          
        
      }
  
VBObox2.prototype.draw = function() {
    //=============================================================================
    // Send commands to GPU to select and render current VBObox contents.
      // check: was WebGL context set to use our VBO & shader program?
      if(this.isReady()==false) {
            console.log('ERROR! before' + this.constructor.name + 
                  '.draw() call you needed to call this.switchToMe()!!');
      }
    
      // ----------------------------Draw the contents of the currently-bound VBO:
      //gl.drawArrays(gl.TRIANGLE_STRIP, sphStart / floatsPerVertex, sphVerts.length / floatsPerVertex);
    }
VBObox2.prototype.reload = function() {
//=============================================================================
// Over-write current values in the GPU for our already-created VBO: use 
// gl.bufferSubData() call to re-transfer some or all of our Float32Array 
// 'vboContents' to our VBO, but without changing any GPU memory allocations.
  											
 gl.bufferSubData(gl.ARRAY_BUFFER, 	// GLenum target(same as 'bindBuffer()')
                  0,                  // byte offset to where data replacement
                                      // begins in the VBO.
 					 				this.vboContents);   // the JS source-data array used to fill VBO
}






function drawParts(gl,ModelMatrix, MvpMatrix,NormalMatrix, u_NormalMatrixLoc,u_MvpMatrixLoc, u_ModelMatrixLoc) {
	//===========================================================

  pushMatrix(ModelMatrix);
  MvpMatrix.set(g_worldMat);
	ModelMatrix.translate(-0.2, 0.8, 0.0);
	//ModelMatrix.rotate(currentAngle, 1, 0, 0);  // spin around y axis.
	ModelMatrix.rotate(90.0, 1.0, 0.0, 0);
	ModelMatrix.scale(0.7, 0.7, 0.7);
  NormalMatrix.setInverseOf(ModelMatrix);
  NormalMatrix.transpose();

  // *: Update our world matrix with ModelMatrix values
  MvpMatrix.multiply(ModelMatrix);
  gl.uniformMatrix4fv(u_ModelMatrixLoc,	false, ModelMatrix.elements);									
  gl.uniformMatrix4fv(u_NormalMatrixLoc, false,NormalMatrix.elements);
  gl.uniformMatrix4fv(u_MvpMatrixLoc, false, MvpMatrix.elements);

	// Draw just the the cylinder's vertices:
	gl.drawArrays(gl.TRIANGLE_STRIP,				// use this drawing primitive, and
								cubeStart/floatsPerVertex, // start at this vertex number, and
		cubeVerts.length / floatsPerVertex);	// draw this many vertices.
	
	// CHILDREN, inherit parent transformations, [ rotated 90 deg, scaled, rotating with global control ]
	// first tree
  pushMatrix(ModelMatrix); 
    MvpMatrix.set(g_worldMat);
		// first segment, rigid
		//ModelMatrix.rotate(-90.0,0.0,0.0,1.0)
		ModelMatrix.scale(0.5, 0.5, 0.5);
		ModelMatrix.translate(-0.5, 0.4, 0.0);
    NormalMatrix.setInverseOf(ModelMatrix);
    NormalMatrix.transpose();
  
    // *: Update our world matrix with ModelMatrix values
    MvpMatrix.multiply(ModelMatrix);
    gl.uniformMatrix4fv(u_ModelMatrixLoc,	false, ModelMatrix.elements);									
    gl.uniformMatrix4fv(u_NormalMatrixLoc, false,NormalMatrix.elements);
    gl.uniformMatrix4fv(u_MvpMatrixLoc, false, MvpMatrix.elements);
  
		// Draw just the the cylinder's vertices:
		gl.drawArrays(gl.TRIANGLE_STRIP,				// use this drawing primitive, and
									cylStart/floatsPerVertex, // start at this vertex number, and
			cylVerts.length / floatsPerVertex);	// draw this many vertices.
	  
		// second segment, can rotate on x
		// second segment, can rotate on x
		ModelMatrix.rotate(angleA, 0.0, 0.0, 1.0);
		ModelMatrix.translate(0.0, 0.5, 0.0);
    NormalMatrix.setInverseOf(ModelMatrix);
    NormalMatrix.transpose();
  
    // *: Update our world matrix with ModelMatrix values
    MvpMatrix.multiply(ModelMatrix);
    gl.uniformMatrix4fv(u_ModelMatrixLoc,	false, ModelMatrix.elements);									
    gl.uniformMatrix4fv(u_NormalMatrixLoc, false,NormalMatrix.elements);
    gl.uniformMatrix4fv(u_MvpMatrixLoc, false, MvpMatrix.elements);
  
    pushMatrix(ModelMatrix)
    MvpMatrix.set(g_worldMat);
    ModelMatrix.translate(0.2, 0.0, 0.0);
    ModelMatrix.rotate(-30.0+angleD, 0.0, 0.0, 1.0)
    
    MvpMatrix.multiply(ModelMatrix);
    gl.uniformMatrix4fv(u_ModelMatrixLoc,	false, ModelMatrix.elements);									
    gl.uniformMatrix4fv(u_NormalMatrixLoc, false,NormalMatrix.elements);
    gl.uniformMatrix4fv(u_MvpMatrixLoc, false, MvpMatrix.elements);
  
    // Draw just the the cylinder's vertices:
    gl.drawArrays(gl.TRIANGLE_STRIP,				// use this drawing primitive, and
                  cylStart/floatsPerVertex, // start at this vertex number, and
      cylVerts.length / floatsPerVertex);	// draw this many vertices.
    
    ModelMatrix.translate(0.0, 0.8, 0.0);
    ModelMatrix.scale(1.0, 1.0, 1.0, 0.0);
    for (j = 0; j < 6; i++, j++) {
      ModelMatrix.rotate(15*j, 0.0, 0.0, 1.0);
      MvpMatrix.multiply(ModelMatrix);
      gl.uniformMatrix4fv(u_ModelMatrixLoc,	false, ModelMatrix.elements);									
      gl.uniformMatrix4fv(u_NormalMatrixLoc, false,NormalMatrix.elements);
      gl.uniformMatrix4fv(u_MvpMatrixLoc, false, MvpMatrix.elements);  
      // Draw just the the cylinder's vertices:
      gl.drawArrays(gl.TRIANGLE_STRIP,				// use this drawing primitive, and
                    cylStart/floatsPerVertex, // start at this vertex number, and
        cylVerts.length / floatsPerVertex);	// draw this many vertices.
      
    }
    ModelMatrix = popMatrix();
	ModelMatrix = popMatrix();                               // reset to isolate it from another tree

}


function makeCone() {
  var errColr = new Float32Array([1.0, 0.2, 0.2]);  // Bright-red trouble color.

  var capVerts = 6.0;  // # of vertices around the bottom cap of the cone
  var topRadius = 0.0;  // top radius is 0 for a cone

  // Create a (global) array to hold all of this cone's vertices;
  coneVerts = new Float32Array((4*capVerts) * floatsPerVertex);

  // Create circle-shaped bottom cap of cone at z=-1.0, radius 1.0
  for (v = 0, j = 0; v < (2 * capVerts) - 1; v++, j += floatsPerVertex) {
    if (v % 2 == 0) {
      coneVerts[j] = Math.cos(Math.PI * v / capVerts);  // x
      coneVerts[j + 1] = Math.sin(Math.PI * v / capVerts);  // y
      coneVerts[j + 2] = -1.0;  // z
      coneVerts[j + 3] = errColr[0];
      coneVerts[j + 4] = errColr[1];
      coneVerts[j + 5] = errColr[2];
    } else {
      coneVerts[j] = 0.0;  // x,y,z,w == 0,0,-1,1; centered on z axis at -1.
      coneVerts[j + 1] = 0.0;
      coneVerts[j + 2] = -1.0;
      coneVerts[j + 3] = errColr[0];
      coneVerts[j + 4] = errColr[1];
      coneVerts[j + 5] = errColr[2];
    }
  }
  var side = 1;
  // Create the cone side walls
  for (v = 0; v < 2 * capVerts; v++, j += floatsPerVertex) {
    if (v % 2 == 0) {
      coneVerts[j] = Math.cos(Math.PI * (v) / capVerts);  // x
      coneVerts[j + 1] = Math.sin(Math.PI * (v) / capVerts);  // y
      coneVerts[j + 2] = -1.0;  // ==z  Cone walls,
    if (side == 1) {
      coneVerts[j + 3] = errColr[0];
      coneVerts[j + 4] = errColr[1];
      coneVerts[j + 5] = errColr[2];
      side=0
    } else {
      coneVerts[j + 3] = errColr[0];
      coneVerts[j + 4] = errColr[1];
      coneVerts[j + 5] = errColr[2];
      side=1
    }

    } else {
      coneVerts[j] = topRadius * Math.cos(Math.PI * (v - 1) / capVerts);  // x
      coneVerts[j + 1] = topRadius * Math.sin(Math.PI * (v - 1) / capVerts);  // y
      coneVerts[j + 2] = 1.0;  // ==z  Cone top,
      coneVerts[j + 3] = errColr[0];
      coneVerts[j + 4] = errColr[1];
      coneVerts[j + 5] = errColr[2];
    }
  }


}
