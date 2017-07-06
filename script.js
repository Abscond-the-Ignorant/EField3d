/*
TODO
Make it so I can add/move charges
Make it is so the node rotates when hit 
	-RotatePhi
Add a prespective
Optimize SeedNodes
*/
var TWO_PI = Math.PI * 2;

var source_lines_per_unit_charge = 6;
var k = 10; // 1/4 pi epsilon naught

// configuration:
var step = 5;
var start_step = 0.001;
var max_steps = 1000;
var Utolerance = 0.001;
var step_equi = 0.1;
var max_equi_step = 100;
var potential_multiple = 3;
var hide_charge_values = false;

var total_charge = 0;
var max_x =  -1e20;
var min_x =  1e20;
var max_y =  -1e20;
var min_y =  1e20;
this.trans = Matrix.I(4);
this.chargeSelected = false;


$(function(){
	gPolar3d = new EField3d($('#viewport'));
	// $('#run_checkbox').change(doTimer);
	/*	
  var self = this;
  $(window).bind('resize',function(ev) { return self.Resize(ev); });
  
  $(window).bind('mousemove',function(ev) { return self.DoMouse(ev); });
  $(this.element).bind('mousedown',function(ev) { return self.DoMouse(ev); });
  $(window).bind('mouseup',function(ev) { return self.DoMouse(ev); });
  $(this.element).bind('mouseout' ,function(ev) { return self.DoMouse(ev); });  
  $('.addcharge').bind('mousedown' ,function(ev) { return self.AddCharge(ev); });  

  $(this.element).bind('touchstart',function(ev) { return self.DoMouse(ev); });
  $(window).bind('touchmove',function(ev) { return self.DoMouse(ev); });
  $(window).bind('touchend',function(ev) { return self.DoMouse(ev); });
  $('.addcharge').bind('touchstart' ,function(ev) { return self.AddCharge(ev); });  

  $('#ctl-do-eqipotential').click(function(){ self.Draw();});
  $('#ctl-zoom-in').click(function(){ self.DoZoom(1); });
  $('#ctl-zoom-out').click(function(){ self.DoZoom(-1); });
		
		
		
	*/
	$('#lines_per_unit_charge').html(source_lines_per_unit_charge);  
	$('#lines_slider').slider({
    	value: source_lines_per_unit_charge,
    	min: 3,
    	max: 26,
    	step: 1,
    	slide: function(event,ui) {  source_lines_per_unit_charge = ui.value; 
    	                            $('#lines_per_unit_charge').html(source_lines_per_unit_charge);
									gPolar3d.StartDraw();
    	                          }
  });
  


  $('#downloadlink').bind('click' ,function(ev) { 
    
    var dt = applet.canvas .toDataURL('image/png');
    this.href = dt;
    
    // return DoPrint($('#everything'),true);  
  });
  
  

  	gPolar3d.charges = [];
  	gPolar3d.StartDraw();
  // doTimer();
});
EField3d.prototype = new Pad3d;           
EField3d.prototype.constructor = EField3d;
function EField3d( element, options ){
  // console.log('TriDView ctor');
  if(!element) {
    // console.log("TriDView: NULL element supplied.");
    return;
  }
  if($(element).length<1) { 
    // console.log()
    return;   
  }
  
  var settings = {
    default_look_at:    [0,0,0],
    default_camera_distance: 800,
    camera_distance_max: 8000,
    camera_distance_min: 50,
    default_theta: -0.1,
    default_phi: 0.5,
  }
  $.extend(true,settings,options);  // Change default settings by provided qualities.
  Pad3d.call(this, element, settings); // Give settings to Pad contructor.
  this.ResetView();
  this.gSetupDirty = true;
  
}
EField3d.prototype.StartDraw = function(){
	this.objects = [];
    // this.charges = [];
    if(this.charges.length==0) this.charges = [
                                { q :  1, x : 150,  y: -50, z: 80},
                                { q : -1, x : 51,  y: 50,  z: -150},
                                { q :  1, x : -151, y: -150, z: 200},
                                { q : -1, x : -50, y: 50,  z: 100},
                                ];



    // if(this.charges.length==0) this.charges = [
    //                             { q :  1, x : 100,  y: -20, z: 20},
    //                             { q : -2, x : -100,   y: 20, z: -20},
    //                             { q : -1, x : 50,   y:-200, z: 170},
    //                             { q :  2, x : -150,   y:80, z: 200},
    //                             ];



	this.DrawCharges();
    this.FindFieldLines();
	this.DrawFieldLines();
	this.Draw();
	
	
}
EField3d.prototype.DrawCharges = function(){
	this.FindTranslationalMatrix();
    for(var i=0 ;i<this.charges.length; i++) {
      var charge = this.charges[i];
	  console.log("q: ",charge.q);
	  this.FindChargePosition(charge);
      total_charge += charge.q;
	  charge.times_seeded = 0;
      charge.r = 10*Math.abs(charge.q);
      charge.n_nodes = Math.round(Math.abs(source_lines_per_unit_charge*charge.q));
      charge.nodes = []; // All successful or unsuccesful nodes
      charge.nodesUsed = []; // Nodes that have actually worked.
      charge.nodesNeeded = []; // Some idea what nodes we should try.
      charge.collisions = []; // where the charge was hit
      if(charge.x > max_x) max_x = charge.x;
      if(charge.x < min_x) min_x = charge.x;
      if(charge.y > max_y) max_y = charge.y;
      if(charge.y < min_y) min_y = charge.y;
	  
	  if(charge.q>0){
		  this.AddPoint(charge.x,charge.y,charge.z,charge.r,"RED","ORANGE",charge);
	  }else if(charge.q<0){
		  this.AddPoint(charge.x,charge.y,charge.z,charge.r,"BLUE","GREEN",charge);
	  }
    }
    this.charges.sort(chargesort);
    if(total_charge<0) this.charges.reverse();
	this.Draw();
}
EField3d.prototype.FindTranslationalMatrix = function(){
    var trans = Matrix.I(4);
    trans = this.Translate(trans,-this.look_at[0],-this.look_at[1],-this.look_at[2]);
    trans = this.RotateY(trans,this.phi);
    trans = this.RotateX(trans,this.theta);
    trans = this.Translate(trans,0,0,this.camera_distance);
	this.trans = trans;
}
EField3d.prototype.FindChargePosition = function(charge){
	//Must call this.FindTranslationalMatrix() first
	var p1 = this.trans.x(Vector.create([charge.x,charge.y,charge.z,1]));
    charge.u =  p1.e(1)/p1.e(3)*this.proj_dist;
    charge.v =  p1.e(2)/p1.e(3)*this.proj_dist;
	charge.tranz = p1.e(3);
}
EField3d.prototype.SetXYCoords = function(charge,u,v){
	console.log("Hi.");
	var XYZVector = Vector.create([u / this.proj_dist * charge.tranz,v / this.proj_dist * charge.tranz,charge.tranz,1]);
	console.log("I'm not bitter.");
	var inverseTrans = this.trans.inverse();
	console.log("You getting here?");
	if(inverseTrans){
		console.log("inverseTrans is not null, ",inverseTrans);
		var a = inverseTrans.x(XYZVector);
		charge.x = a.e(1);
		charge.y = a.e(2);
		charge.z = a.e(3);
		this.FindChargePosition(charge);
	}
	
}
EField3d.prototype.FindFieldLines = function(){
	
    this.fieldLines = [];
	
    for(var i=0 ;i<this.charges.length; i++) {
      var charge = this.charges[i];
      // console.log("Find field lines for charge ",i," with charge ",charge.q);
      // console.log("Doing charge",i,"with q=",charge.q,"which has ",charge.nodesUsed.length,"/",charge.n_nodes," nodes");


      while(charge.nodesUsed.length < charge.n_nodes && charge.nodes.length<source_lines_per_unit_charge*5) {
        if(charge.nodes.length>source_lines_per_unit_charge*4) {
          console.warn("Wow! Tried way too many nodes.",charge.nodes);
        }
      // console.log("Doing charge",i,"with q=",charge.q,"which has ",charge.nodesUsed.length,"/",charge.n_nodes," nodes");
        // console.log("Doing node on charge",i);
      
        var start_angle = this.FindNodePosition(charge);
        var r = charge.r;
        // Boost in initial direction by radius.
        var fieldline = { startCharge: charge };
        fieldline.start = "charge";
		
		
		fieldline.start_x = charge.x + charge.r*Math.sin(start_angle.phi)*Math.cos(start_angle.theta);
		fieldline.start_y = charge.y + charge.r*Math.sin(start_angle.phi)*Math.sin(start_angle.theta);
		fieldline.start_z = charge.z + charge.r*Math.cos(start_angle.phi);
        fieldline.start_angle = start_angle;
		
        var dir = 1;
        if(charge.q<0) dir = -1;
        fieldline.dir     = dir;
      
        var nodeFinished = this.TraceFieldLine(fieldline); 
        if(nodeFinished) {
          charge.nodesUsed.push(start_angle);
          this.fieldLines.push(fieldline);      
        }
      } // nodeFinished
    }
}
EField3d.prototype.TraceFieldLine = function(fieldline){
  // console.log(fieldline);
  var x = fieldline.start_x;
  var y = fieldline.start_y;
  var z = fieldline.start_z;
  
  fieldline.points  = [{x:x,y:y,z:z}];
  var lastE = this.Field(x,y,z);

  var traceFinished = false;
  var nstep = 0;
  while(true) {
    nstep++;
    var E = this.Field(x,y,z);
    
    var dx = E.gx * step *fieldline.dir;
    var dy = E.gy * step *fieldline.dir;
    var dz = E.gz * step *fieldline.dir;
    
    // This check doesn't work.
    // I've also tried the form of E/U, which should check for small field in high-potential areas (which is what I want to check)
    // I can't just check for small field, because that happens naturally at large distances from the middle.
    // if(Math.abs(E.E) < 1 && Math.abs(E.U) > 2) {
    //   console.log("Line went through effective zero-field region. This isn't a good line. E=",E.E," U=",E.U, " step=",step);
    //   return false;
    // }
    
    x += dx;
    y += dy;
    z += dz;
    fieldline.points.push({x:x,y:y,z:z});
    lastE = E;
    var collide = this.FindCollision(x,y,z);
    if(collide && (fieldline.dir*collide.q < 0) && nstep>1) {
      // Find the best possible node for this line.
      if(collide.n_nodes >= collide.nodes.length+1 == 0) {
        // console.log("Line failed - hit q=",collide.q,"which has no nodes left.");
        // Comment these lines out if you want it to just sail through without stopping...
        // console.warn("Line failed - hit q=",collide.q,"which has no nodes left.");
        return false; //nodeFinished=false;
      } else {
        this.DoCollision(collide,x,y,z);
        fieldline.endCharge = collide;
        fieldline.nstep = nstep;
        // console.log("Line succeeded - hit q=",collide.q);
        return true; // nodeFinished
      }
    }
	
    if(nstep>max_steps){
      fieldline.endCharge = null;
      fieldline.endAngle     = null;
      fieldline.endNodeAngle = null;
      fieldline.nstep = nstep;
      // console.log("Line succeeded - no hit");
      return true;
    }  // if nstep
	 
  } // trace loop
}
EField3d.prototype.DrawFieldLines = function(){
	for(var k = 0; k < this.fieldLines.length; k++){
		for(var q = 0; q < this.fieldLines[k].points.length - 1; q++){
			var p = this.fieldLines[k].points;
			// console.log(this.fieldLines[k]);
			this.AddLine(p[q].x,p[q].y,p[q].z,p[q+1].x,p[q+1].y,p[q+1].z,1,"BLACK",null);
		}
		this.AddArrows(this.fieldLines[k]);
	}
}
EField3d.prototype.FindNodePosition = function(charge){
  // If this is a virgin charge. Find seed positions.
  if(charge.nodes.length == 0 && charge.nodesNeeded.length == 0) {
    this.SeedNodes(charge,{theta: 0, phi: 0});
  }

  // See if there are some precomputed positions to try.
  if(charge.nodesNeeded && charge.nodesNeeded.length>0) {
    var t = charge.nodesNeeded.shift();
	// console.log(t);
    charge.nodes.push(t);
    return t;
  }
  
//try to seed it again,this time with more nodes
  this.SeedNodes(charge,{theta: 0, phi: 0});
  if(charge.nodesNeeded && charge.nodesNeeded.length>0) {
    var t = charge.nodesNeeded.shift();
    charge.nodes.push(t);
    return t;
  }
  
}
EField3d.prototype.SeedNodes = function(charge,startangle){
	var n = 0;
	var q = 0;
	while(charge.n_nodes > 2 * (n*n - n + 1)){ n++; }
	n += charge.times_seeded;
	charge.times_seeded++;
	for(var i = 0; i<2 * (n*n - n + 1); i++) {
	  if(((i+q) % (1 * n) == 0) && i > n){ q++; }
	  var theta = Math.floor((i+q)/(2*n)) * Math.PI/n;
  	  var phi = ((i + q) %(2*n)) * Math.PI/n;
	  // console.log("phi: ",phi);
	  // console.log("theta: ",theta);
	  
	  /* 
	  var a = {x: r*Math.sin(phi)*Math.cos(theta) ,y: r*Math.sin(phi)*Math.sin(theta) ,z: r*Math.cos(phi)};
	  var b = RotateTheta(
		  {x: r*Math.sin(phi)*Math.cos(theta) ,y: r*Math.sin(phi)*Math.sin(theta) ,z: r*Math.cos(phi)},
		  startangle.theta);
	  
	  
	  var a = RotatePhi(RotateTheta(
		  		  {x: r*Math.sin(phi)*Math.cos(theta) ,y: r*Math.sin(phi)*Math.sin(theta) ,z: r*Math.cos(phi)},
			  startangle.theta),startangle.phi);
	  */
	  charge.nodesNeeded.push({theta:theta,phi:phi});
	  // charge.nodesNeeded.push(GetAngle({x: r*Math.sin(phi)*Math.cos(theta) ,y: r*Math.sin(phi)*Math.sin(theta) ,z: r*Math.cos(phi)}));
/*
	  charge.nodesNeeded.push( 
		  GetAngle(
			  RotatePhi(
				  RotateTheta(
					  {x: r*Math.sin(phi)*Math.cos(theta) ,y: r*Math.sin(phi)*Math.sin(theta) ,z: r*Math.cos(phi)},
				  startangle.theta)
			  ,startangle.phi)));
	  */
  }
}
EField3d.prototype.AddArrows = function(fieldline){
	var n = fieldline.points.length;
	if(n > 500){ for(var i = 1; i < Math.floor(n / 200); i++) {if(n - Math.round(n / Math.floor(n / 200)) * i > 50){
		this.AddArrow(fieldline.points[Math.round(n / Math.floor(n / 200)) * i],
		fieldline.points[Math.round(n / Math.floor(n / 200)) * i - fieldline.dir]);
	}}}
	else{
		this.AddArrow(fieldline.points[Math.round(n / 2)], fieldline.points[Math.round(n/2)-fieldline.dir]);
	}
}
EField3d.prototype.AddArrow = function(point,prevPoint){
	// console.log(point);
			// console.log("point: ", point);
			// console.log("prevPoint: ", prevPoint);
	var arrow = Vector.create([point.x - prevPoint.x, point.y - prevPoint.y, point.z - prevPoint.z]);
	if(!arrow.isParallelTo(Vector.i)){
		var p1 = arrow.cross(Vector.i);
		p1 = p1.toUnitVector();
		p1 = p1.x(0.5);
		var arrowLength = 20/step;
		var arrowAxis = Line.create([0,0,0],[arrow.e(1),arrow.e(2),arrow.e(3)]);
		for(var i = 0; i < 4; i ++){
			var arrowLine = Vector.create([prevPoint.x+p1.e(1)-point.x,prevPoint.y+p1.e(2)-point.y,prevPoint.z+p1.e(3)-point.z]);
			this.AddLine(point.x,point.y,point.z,point.x+arrowLength*arrowLine.e(1),point.y+arrowLength*arrowLine.e(2),point.z+arrowLength*arrowLine.e(3),1,"BLACK",null);
			p1 = p1.rotate(2*Math.PI/4,arrowAxis);
		}
	}else if(!arrow.isParallelTo(vector.j)){
		var p1 = arrow.cross(Vector.j);
		p1 = p1.toUnitVector();
		p1 = p1.x(0.5);
		var arrowLength = 5;
		var arrowAxis = Line.create([0,0,0],[arrow.e(1),arrow.e(2),arrow.e(3)]);
		for(var i = 0; i < 4; i ++){
			var arrowLine = Vector.create([prevPoint.x+p1.e(1)-point.x,prevPoint.y+p1.e(2)-point.y,prevPoint.z+p1.e(3)-point.z]);
			this.AddLine(point.x,point.y,point.z,point.x+arrowLength*arrowLine.e(1),point.y+arrowLength*arrowLine.e(2),point.z+arrowLength*arrowLine.e(3),1,"BLACK",null);
			console.log(p1);
			p1 = p1.rotate(2*Math.PI/4,arrowAxis);
		}
	}
}
EField3d.prototype.Field = function(x,y,z){
  var Ex = 0;
  var Ey = 0;
  var Ez = 0;
  var U  = 0;
  var dUdx = 0;
  for(var i=0 ;i<this.charges.length; i++) {
    var c = this.charges[i];
    var dx = x-c.x;
    var dy = y-c.y;
    var dz = z-c.z;
    var r2 = dx*dx+dy*dy+dz*dz;
    var r = Math.sqrt(r2);
    var E = c.q/r2;
    Ex += dx/r*E;
    Ey += dy/r*E;
    Ez += dz/r*E;
    U += c.q/r;
  }
  var E2 = Ex*Ex + Ey*Ey + Ez*Ez;
  var E = Math.sqrt(E2);
  var ret = { x: x, y: y, z: z,  			// Coordinates.
              U: U,              			// Potential
              E: E,              			// Field magnitude
	  		  Ex: Ex, Ey: Ey, Ez: Ez,    	// Field vector
		  	  gx: Ex/E, gy: Ey/E, gz: Ez/E 	// Field direction vector
              };
  // console.log("Field at "+x+","+y,ret);
  return ret;
}
EField3d.prototype.FindCollision = function(x,y,z){
  for(var i=0 ;i<this.charges.length; i++) {
    var c = this.charges[i];
    var dx = x-c.x;
    var dy = y-c.y;
    var dz = z-c.z;
    var r2 = dx*dx+dy*dy+dz*dz;
    var cr = c.r-0.01;
	// console.log(c.x,c.y,c.z);
	// console.log(x,y,z);
    if(r2 < (cr*cr)) {
      // console.log("collision",dx,dy,r2,cr*cr);
      return c;
    }
  }
  return null;
}
EField3d.prototype.DoCollision = function(collide,x,y,z){
  // console.warn("collided with charge that has ",collide.nodesNeeded.length,"left ")
  collide.collisions.push({x:x,y:y,z:z});
  dx = x-collide.x;
  dy = y-collide.y;
  dz = z-collide.z;
  var angle = GetAngle({x:dx,y:dy,z:dz});
  
  
  collide.nodes.push(angle);
  collide.nodesUsed.push(angle);
  
  if(collide.nodesUsed.length==1) {
    // This is the first line to collide. Seed other positions around this.
	  
	  /*############
	  ##############
	  ##############
	  ############*/
    this.SeedNodes(collide,{phi:0,theta:0});
    // this.SeedNodes(collide,angle);
  }
    
  var best = 0;
  var bestdiff = 9e9;
  for(var k=0; k<collide.nodesNeeded.length;k++){
      var diff = GetDiff(collide.nodesNeeded[k],angle);
      if(diff<bestdiff) {bestdiff = diff; best = k};
  }
  collide.nodesNeeded.splice(best,1);
}
function RotateTheta(a, angle){
	var b = {
		x: a.x*Math.cos(angle)-a.y*Math.sin(angle),
		y: a.x*Math.sin(angle)-a.y*Math.cos(angle),
		z: a.z
	};
	return b;
}
function RotatePhi(a, angle){
	if(angle == 0) return a;
	// console.log("x: ", a.x);
	// console.log("y: ", a.y);
	var k = {x:-a.y / Math.sqrt(a.x*a.x + a.y*a.y),y: a.x / Math.sqrt(a.x*a.x + a.y*a.y), z: 0};
	var dot = k.x*a.x + k.y*a.y;
	var c = Math.cos(angle);
	var s = Math.sin(angle);
	var b = {
		x: c*a.x + dot*(1-c)*k.x + s*k.y*a.z,
		y: c*a.y + dot*(1-c)*k.y - s*k.z*a.z,
		z: c*a.z + (k.x*a.y-k.y*a.x)*s
	};
	return b;
}
function chargesort(a,b){
  var cmp = a.q - b.q;
  if(cmp==0) cmp = a.y - b.y;
  return cmp;
}
function GetDiff(angle1,angle2){
	var a = {x: Math.sin(angle1.phi)*Math.cos(angle1.theta), y: Math.sin(angle1.phi)*Math.sin(angle1.theta), z: Math.cos(angle1.phi)};
	var b = {x: Math.sin(angle2.phi)*Math.cos(angle2.theta), y: Math.sin(angle2.phi)*Math.sin(angle2.theta), z: Math.cos(angle2.phi)};
	return Math.abs( Math.acos( Dot(a,b) ) % (TWO_PI) );
}
function GetAngle(a){ return {theta: Math.atan(a.y/a.x), phi: Math.atan(Math.sqrt(a.x*a.x+a.y*a.y)/a.z)}; }	
function Dot(a,b){ return a.x * b.x + a.y * b.y + a.z * b.z; }
EField3d.prototype.mouseMove = function(ev){
	// if(this.dragging){
	// this.drag(ev);
	//   	// return;
	//   }
  // Called for any movement in pad.  
  // find highlighed object.
  // go through object list from back to front, so that the frontmost object will highlight.
  var offset = getAbsolutePosition(this.canvas);
  var x = ev.pageX - offset.x;
  var y = ev.pageY - offset.y;    
  var u = this.width/2 -x;  // Translate and flip inverse to screen coords
  var v = this.height/2-y;

  var selected = null;
  for(var i=0;i<this.objects.length;i++)
  {
    var obj = this.objects[i];
    // Only compute if object is selectable: i.e. it has a reference source (is not scenery)
    if(obj.source){
      if(obj.type=='l') {
        // if(sample++%100==0) console.log(obj.source,
        //                                 u,v,obj.au,obj.av,obj.bu,obj.bv,
        //                                 GeoUtils.line_to_point(u,v,obj.au,obj.av,obj.bu,obj.bv));
        if( GeoUtils.line_is_close_to_point(u,v,obj.au,obj.av,obj.bu,obj.bv,
          this.mouse_highlight_range*this.proj_dist/obj.meanz) )
          selected = obj.source; 
        }
      }
      if(obj.type=='p') {
		  if(obj.selected){
	  		obj.selected = false;
			this.Draw();
		  }
        var du = u - obj.au;
        var dv = v - obj.av;
        var d = Math.sqrt(du*du+dv*dv);
        if( d < this.mouse_highlight_range*this.proj_dist/obj.meanz ){
			obj.selected = true;
			this.DrawObj(obj);
        	selected = obj.source;
        }
      }
  }
  if(selected) this.HoverObject(selected);
	if(this.dragging){
	this.drag(ev);
  	// return;
  }
}
EField3d.prototype.drag = function(ev){
	  // Called for any movement in page.
	  if(!this.dragging)  return;  
	   var x = ev.pageX;
	   var y = ev.pageY;
	   console.log("dragging....."+ev.type+" "+x+" "+y);
	   if(ev.originalEvent.touches) {
	     if(ev.originalEvent.touches) console.log("Touches avalable");
	     if (ev.originalEvent.touches.length > 1) {this.dragging = false; return; } // don't allow multi-touch
	     x = ev.originalEvent.touches[0].pageX;
	     y = ev.originalEvent.touches[0].pageY;       
	     ev.originalEvent.preventDefault();       
	     console.log("Touch: " + x + " " + y)
	   } else {
	     if(this.startDragTouch) return; // Avoid processing a fake mousemove event when running ios.
	   }
	   if(!this.chargeSelected){
	   		var dx = (x - this.startDragX);
	   		var dy = (y - this.startDragY);
	   		console.log("Rotate: " + dx+" "+dy);
	   		
	   		var sx =  dx / $(this.element).width();
	   		var sy =  dy / $(this.element).height();
	   		console.log("Rotate: " + sx+" "+sy);
	   		
	   		this.startDragX = x;
	   		this.startDragY = y;
	   		
	   		if((this.mouse_mode == "pan" && this.drag_button != 3)||(this.mouse_mode == "rotate" && this.drag_button == 3)) {
	   		  // panning
	   		
	   		  // Need to find 'up' and 'right' vectors to find how to shift look-at point.
	   		  var trans = Matrix.I(4);
	   		  trans = this.RotateX(trans,-this.theta);
	   		  trans = this.RotateY(trans,-this.phi);
	   		  var u = dx;
	   		  var v = dy;
	   		  var inv = Vector.create([u,v,0,1]);
	   		  var outv = trans.x(inv);
	   		  this.panMe(null,outv.e(1),outv.e(2),outv.e(3));
	   		  
	   		} else {
	   		  // rotating
	   		  sx = Math.asin(sx);
	   		  sy = -Math.asin(sy);
	   		  if(!ev.originalEvent.touches) {
	   		    sx=5*sx;
	   		    sy=5*sy;
	   		  }
	   		  console.log("Rotate: " + sx+" "+sy);
	   		  this.rotateMe(sy,sx);
	   		}
	   		this.Draw();
   	}else{
		var offset = getAbsolutePosition(this.canvas);
		var x = ev.pageX - offset.x;
		var y = ev.pageY - offset.y;    
		var u = this.width/2 -x;  // Translate and flip inverse to screen coords
		var v = this.height/2-y;
		this.FindTranslationalMatrix();
	    for(var i=0;i<this.charges.length;i++){
	    	var charge = this.charges[i];
			this.FindChargePosition(charge);
	        var du = u - charge.u;
	        var dv = v - charge.v;
	        var d = Math.sqrt(du*du+dv*dv);
	        if( d < charge.r*this.proj_dist/charge.tranz*1.5 ){
				this.SetXYCoords(charge,u,v);
				this.StartDraw();
	        }
	    }
		
   	}
}
EField3d.prototype.startDragging = function(ev){
  for(var i=0;i<this.objects.length;i++) if(this.objects[i].selected) this.chargeSelected = true;
  
  // this.startDragTransform = this.transform_matrix.dup();
  this.startDragX = ev.pageX;
  this.startDragY = ev.pageY;
  if(ev.originalEvent.touches) {
    if (ev.originalEvent.touches.length > 1) {this.dragging = false;return; } // don't allow multi-touch
    this.startDragX = ev.originalEvent.touches[0].pageX;
    this.startDragY = ev.originalEvent.touches[0].pageY;       
    this.startDragTouch = true;
    ev.originalEvent.preventDefault();       
  }
  /*
  	for(var i=0;i<this.objects.length;i++){
	  var obj = this.objects[i];
	  if(obj.selected){


  		var offset = getAbsolutePosition(this.canvas);
  		var x = ev.pageX - offset.x;
  		var y = ev.pageY - offset.y;
  		var u = this.width/2 -x;  // Translate and flip inverse to screen coords
  		var v = this.height/2-y;
  		this.FindTranslationalMatrix();
  	    for(var i=0;i<this.charges.length;i++){
  	    	var charge = this.charges[i];
  			this.FindChargePosition(charge);
  	        var du = u - charge.u;
  	        var dv = v - charge.v;
  	        var d = Math.sqrt(du*du+dv*dv);
  	        if( d < charge.r*this.proj_dist/charge.tranz*1.5 ){
  				this.SetXYCoords(charge,u,v);
  				this.StartDraw();
  	        }
  	    }
		
		
	  }
	}*/
  console.log("start drag "+this.startDragX+" "+this.startDragY);
  this.dragging = true;
  this.drag_button = ev.which;
  $(this.element).css("cursor","move !important");
  return true;
  
}
EField3d.prototype.stopDragging = function(ev){
  $(this.element).css("cursor","auto");
  this.chargeSelected = false;
  this.dragging = false;
}
/*
	code for doing escaping lines -- was in FindFieldLines --


    var escaping_lines = Math.abs(total_charge* source_lines_per_unit_charge);
	
	var n = 1;
	var q = 0;
	while(escaping_lines > 2 * (n*n - n + 1)){ n++; }
      var r = 10000;
    for(var i=0;i<escaping_lines;i++) {
      console.log("Doing escaping line.");
      // Find a position very far away from the charges.
	  if(((i+q)%(n) == 0) && i > n){ q++; }
	  // if(((i+q)%(n*2) == 0) && Math.floor((i+q)/(2*n)) != 0){ q++; }
	  var theta = Math.floor((i+q)/(2*n)) * Math.PI/n;
	  var phi = ((i + q) %(2*n)) * Math.PI/n;

	  var x =  r*Math.sin(phi)*Math.cos(theta);
	  var y =  r*Math.sin(phi)*Math.sin(theta);
 	  var z =  r*Math.cos(phi);
 	  var fieldline = { startCharge: null };
	  fieldline.start_x = x;
	  fieldline.start_y = y;
	  fieldline.start_z = z;
	  if(total_charge > 0)  fieldline.dir = -1;
	  else                  fieldline.dir = 1;
	  fieldline.start = "outside";
	  var nodeFinished = this.TraceFieldLine(fieldline);
	  if(nodeFinished) {
	    this.fieldLines.push(fieldline);      
	  } else {
	    console.log("incoming line failed");
	  }
    
    }
	
*/
/*
function Applet(element, options){
  if(!element) { 
    console.log("Pad: NULL element provided."); return; 
  }
  if($(element).length<1) { 
    console.log("Pad: Zero-length jquery selector provided."); return;
  }
  this.element = $(element).get(0); 
  
  this.bg_color = "255,255,255";
  this.origin_x = 0.0;
  this.origin_y = 0.0;
  this.width_x  = 10.0;
  this.dragging = false;
  
  // Merge in the options.
  $.extend(true,this,options);

  // Merge in options from element
  var element_settings = $(element).attr('settings');
  var element_settings_obj={};
  // if(element_settings) {
  //   eval( "var element_settings_obj = { " + element_settings + '};');; // override from 'settings' attribute of html object.
  //   console.log(element_settings, element_settings_obj);
  //   $.extend(true,this,element_settings_obj); // Change default settings by provided overrides.
  //
  // }
  // Look for an existing canvas, and build one if it's not there.
  if($('canvas',this.element).length<1) {
    this.canvas = document.createElement("canvas");
    this.element.appendChild(this.canvas);
  } else {
    this.canvas = $('canvas',this.element).get(0);    
  }
  
  if(!element) { 
    console.log("Pad: NULL element provided."); return; 
  }
  if($(element).length<1) { 
    console.log("Pad: Zero-length jquery selector provided."); return;
  }
  this.element = $(element).get(0); 


  // Build the drawing context.
  this.ctx = this.canvas.getContext('2d');
  // if(initCanvas) this.ctx = initCanvas(this.canvas).getContext('2d');
  if(!this.ctx) console.log("Problem getting context!");
  if( !$(this.element).is(":hidden") ) {
    width = $(this.element).width();
    height = $(this.element).height(); 
  }
  this.canvas.width = this.width = width;
  this.canvas.height = this.height = height;

  

  // Data.
  this.charges = [];
  
  // see if there's an override in the URL
  var urlParams;
  var match,
      pl     = /\+/g,  // Regex for replacing addition symbol with a space
      search = /([^&=]+)=?([^&]*)/g,
      decode = function (s) { return decodeURIComponent(s.replace(pl, " ")); },
      query  = window.location.search.substring(1);

  urlParams = {};
  while (match = search.exec(query))
         urlParams[decode(match[1])] = decode(match[2]);
  
  for( p in urlParams ) {
    if(p.match(/^q/)) {
      var list = urlParams[p].split(',');
      this.charges.push({q:parseFloat(list[0]),
                         x:parseFloat(list[1]),
                         y:parseFloat(list[2]),
                         r:Math.abs(parseFloat(list[0]))*0.12});
    }
  }
  if(urlParams.lines) {
    source_lines_per_unit_charge = urlParams.lines;
  }
  if(urlParams.hideq) {
    hide_charge_values = true;
  }
  
  if(this.charges.length==0) this.charges = [
                              { q : 1,  x : -1,  y: 1 , r:0.12},
                              { q : -1, x : 1,   y:-0 , r:0.12 },
                              { q : -2, x : 1.001,   y:-1 , r:Math.sqrt(2)*0.12 },
                              ];
  
  
  this.FindFieldLines();
  this.Draw();

  var self = this;
  $(window).bind('resize',function(ev) { return self.Resize(ev); });
  
  $(window).bind('mousemove',function(ev) { return self.DoMouse(ev); });
  $(this.element).bind('mousedown',function(ev) { return self.DoMouse(ev); });
  $(window).bind('mouseup',function(ev) { return self.DoMouse(ev); });
  $(this.element).bind('mouseout' ,function(ev) { return self.DoMouse(ev); });  
  $('.addcharge').bind('mousedown' ,function(ev) { return self.AddCharge(ev); });  

  $(this.element).bind('touchstart',function(ev) { return self.DoMouse(ev); });
  $(window).bind('touchmove',function(ev) { return self.DoMouse(ev); });
  $(window).bind('touchend',function(ev) { return self.DoMouse(ev); });
  $('.addcharge').bind('touchstart' ,function(ev) { return self.AddCharge(ev); });  

  $('#ctl-do-eqipotential').click(function(){ self.Draw();});
  $('#ctl-zoom-in').click(function(){ self.DoZoom(1); });
  $('#ctl-zoom-out').click(function(){ self.DoZoom(-1); });
  

}
Applet.prototype.DoZoom = function( zoom ){
  this.width_x -= zoom;
  this.Draw();
}
Applet.prototype.Resize = function(){
  console.log("Applet::Resize()",this);
  var width = $(this.element).width();
  var height = $(this.element).height(); 
  this.canvas.width = this.width = width;
  this.canvas.height = this.height = height;
  this.Draw();
}
Applet.prototype.Clear = function(){
  //console.log("Pad.Clear()");
  if (!this.ctx) return;
  this.ctx.fillStyle = "rgb("+this.bg_color+")";
  this.ctx.fillRect(0,0,this.canvas.width,this.canvas.height);
}
Applet.prototype.Field = function(x,y){
  var Ex = 0;
  var Ey = 0;
  var U  = 0;
  var dUdx = 0;
  for(var i=0 ;i<this.charges.length; i++) {
    var c = this.charges[i];
    var dx = x-c.x;
    var dy = y-c.y;
    var r2 = dx*dx+dy*dy;
    var r = Math.sqrt(r2);
    var E = c.q/r2;
    Ex += dx/r*E;
    Ey += dy/r*E;
    U += c.q/r;
  }
  var E2 = Ex*Ex + Ey*Ey;
  var E = Math.sqrt(E2);
  var ret = { x: x, y: y,         // Coordinates.
              U: U,              // Potential
              E: E,              // Field magnitude
              Ex: Ex, Ey: Ey,    // Field vector
              gx: Ex/E, gy: Ey/E // Field direction vector
              };
  // console.log("Field at "+x+","+y,ret);
  return ret;
}
Applet.prototype.FindCollision = function(x,y){
  for(var i=0 ;i<this.charges.length; i++) {
    var c = this.charges[i];
    var dx = x-c.x;
    var dy = y-c.y;
    var r2 = dx*dx+dy*dy;
    var cr = c.r-0.0001;
    if(r2 < (cr*cr)) {
      // console.log("collision",dx,dy,r2,cr*cr);
      return c;
    }
  }
  return null;
}
function SpansIntegerMultiple(a,b,r){
  // Does (a,b) span an a value that is an integer multiple of r?
  var da = Math.floor(a/r);
  var db = Math.floor(b/r);
  if(da==db) return null;
  return Math.max(da,db);

}
function PointTripletOrientation(p,q,r){
  // From http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    // See 10th slides from following link for derivation of the formula
    // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    var val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
 
    if (val == 0) return 0;  // colinear
 
    return (val > 0)? 1: 2; // clock or counterclock wise
}
function PointOnSegment( p, q, r){
    if (q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x) &&
        q.y <= Math.max(p.y, r.y) && q.y >= Math.min(p.y, r.y))
       return true;
 
    return false;
}
function LineSegmentsIntersect(p1,q1,  // first line segment points
                               p2,q2   // second line segment points
                            ){
  // From http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
  
  // Find the four orientations needed for general and
  // special cases
  var o1 = PointTripletOrientation(p1, q1, p2);
  var o2 = PointTripletOrientation(p1, q1, q2);
  var o3 = PointTripletOrientation(p1, p2, q2);
  var o4 = PointTripletOrientation(q1, p2, q2);

  var d2 = (q2.x-p1.x)*(q2.x-p1.x) + (q2.y-p1.y)*(q2.y-p1.y)
  // console.log("Check for intersection",o1,o2,o3,o4,Math.sqrt(d2),p1,q1,p2,q2);
  // General case
  if ((o1 != o2) && (o3 != o4))
      return true;

  // Not important to us; tested segments should always be nearly perpendicular
  // Special Cases: check for colinearity.
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  // if (o1 == 0 && PointOnSegment(p1, p2, q1)) return true;
  //
  // // p1, q1 and p2 are colinear and q2 lies on segment p1q1
  // if (o2 == 0 && PointOnSegment(p1, q2, q1)) return true;
  //
  // // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  // if (o3 == 0 && PointOnSegment(p2, p1, q2)) return true;
  //
  //  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  // if (o4 == 0 && PointOnSegment(p2, q1, q2)) return true;

  return false; // Doesn't fall in any of the above cases
}
Applet.prototype.FindNodePosition = function(charge){
  // If this is a virgin charge. Find seed positions.
  if(charge.nodes.length == 0 && charge.nodesNeeded.length == 0) { 
    this.SeedNodes(charge,0);
  }

  // See if there are some precomputed positions to try.
  if(charge.nodesNeeded && charge.nodesNeeded.length>0) {
    var t = charge.nodesNeeded.shift();
    charge.nodes.push(t);
    return t;
  }


  charge.nodes.sort();
  // bifurcate biggest arc
  var biggest_gap = 0;
  var gap_after = 0;
  for(var i=0;i<charge.nodes.length;i++) {
    var t1 = charge.nodes[i];
    var t2;
    if(i+1 < charge.nodes.length) t2 = charge.nodes[(i+1)];
    else t2 = charge.nodes[(i+1)%charge.nodes.length]+TWO_PI; // wrap around
    var dt = Math.abs(t2-t1);
    // console.log(i,t1,t2,dt);
    if(dt>biggest_gap) { gap_after = i; biggest_gap = dt; }
  }
  var new_node = (charge.nodes[gap_after] + biggest_gap*0.5)%(TWO_PI);
  charge.nodes.push(new_node);
  return new_node;
}
Applet.prototype.FindPositionOfU  = function(input,Utarget,Utolerance){
  // Takes an input object: {E: E{E,U,x,y}, x, y}}
  // Returns a similar output object at the best guess for Utarget.
  // Follows the field line at point input.x,input.y until it converges on the Utarget
  
  // We know that U = - delE
  // so by newton-raphson method...
  // distance to go along field line = (U1 - Utarget) / E
  var out = input;
  var it = 0;
  while(Math.abs(out.U - Utarget) > Utolerance) {
    it++;
    var delta = (out.U - Utarget) / out.E;
    var x = out.x + ( delta * out.gx );
    var y = out.y + ( delta * out.gy );
    out = this.Field(x,y);
  }
  // console.log("converge in ", it, "Accuracy: ",out.E.U - Utarget);
  return out;
}
Applet.prototype.SeedNodes = function(charge,startangle){
  // // Original algorithm: Space 'needed' nodes around evenly.
  for(var j = 0; j<charge.n_nodes; j++) {
    charge.nodesNeeded.push(
      (startangle + TWO_PI*j/charge.n_nodes)%TWO_PI 
    );
  }


  // // Algorithm 2: Space 'needed' nodes around accoring to the
  // // LOCAL field, as adjusted by other local charges!
  // var nGrid = 72;
  // var biggestField = 0;
  // var biggestFieldJ = 0;
  // var totField = 0;
  // var grid = [];
  // for(var j=0;j<nGrid;j++) {
  //   var theta = 2*Math.PI*j/nGrid;
  //   var x = charge.x+charge.r*Math.cos(theta);
  //   var y = charge.y+charge.r*Math.sin(theta);
  //   var E = this.Field(x,y);
  //   // console.log(x,y,E,charge);
  //   if(Math.abs(E.E)>biggestField) { biggestField = Math.abs(E.E); biggestFieldJ=j;}
  //   totField += Math.abs(E.E);
  //   grid.push(E.E);
  // }
  // // Now, evenly space them around in integrated field units.
  // var spacing = totField/charge.n_nodes;
  // charge.nodesNeeded.push(2*Math.PI*biggestFieldJ/nGrid);
  //
  // var sum = 0;
  // for(var j=1;j<nGrid;j++) {
  //   var jj = (j+biggestFieldJ)%nGrid;
  //   sum += grid[jj];
  //   if(sum>spacing) {charge.nodesNeeded.push(2*Math.PI*jj/nGrid); sum -= spacing;}
  // }
  // var spacings = [];
  // for(var j=1;j<charge.nodesNeeded.length;j++) {
  //   spacings.push((charge.nodesNeeded[j]-charge.nodesNeeded[j-1])/2/Math.PI);
  // }
  // // console.log('nodes',charge.nodesNeeded);
  // // console.log('spacings',spacings);
  // if(charge.nodesNeeded.length != charge.n_nodes) console.log("Got wrong number of needed points. Wanted ",charge.n_nodes," got ",charge.nodesNeeded.length);

  // Algorithm 3: track from the very center, using epsilon push away from charge center.
  // for(var j = 0; j<charge.n_nodes; j++) {
  //   var theta = 2*j/charge.n_nodes*Math.PI;
  //   var dir = 1;
  //   if(charge.q<0) dir = -1;
  //   var x = charge.x + 0.01*Math.cos(theta);
  //   var y = charge.y + 0.01*Math.sin(theta);
  //   var deltax = 0;
  //   var deltay = 0;
  //   var d2 = 0;
  //   var nstart = 0;
  //   do {
  //     var E = this.Field(x,y);
  //     var dx = E.gx * step/10;
  //     var dy = E.gy * step/10;
  //     x += dx*dir;
  //     y += dy*dir;
  //     deltax = x-charge.x;
  //     deltay = y-charge.y;
  //     d2 = deltax*deltax + deltay*deltay;
  //     nstart++;
  //   } while ( d2 < charge.r*charge.r );
  //
  //   var angle = Math.atan2(deltay,deltax);
  //
  //   charge.nodesNeeded.push(angle);
  //   console.log("need node:",deltay,deltax,angle,nstart);
  //  }

  // charge.nodesRequested = charge.nodesNeeded.slice(0);

}
Applet.prototype.DoCollision = function(collide,x,y){
  // console.warn("collided with charge that has ",collide.nodesNeeded.length,"left ")
  dx = x-collide.x;
  dy = y-collide.y;
  var angle = (Math.atan2(dy,dx)+TWO_PI)%TWO_PI;
  
  
  collide.nodes.push(angle);
  collide.nodesUsed.push(angle);
  
  if(collide.nodesUsed.length==1) {
    // This is the first line to collide. Seed other positions around this.
    this.SeedNodes(collide,angle);
  }
    
  var best = 0;
  var bestdiff = 9e9;
  for(var k=0; k<collide.nodesNeeded.length;k++){
      var diff = Math.abs( (collide.nodesNeeded[k] - angle)%(2*Math.PI) );
      if(diff<bestdiff) {bestdiff = diff; best = k};
  }
  collide.nodesNeeded.splice(best,1);
  
}
Applet.prototype.TraceFieldLine = function(fieldline){
  console.log(fieldline);
  var x = fieldline.start_x;
  var y = fieldline.start_y;
  
  fieldline.points  = [{x:x,y:y}];
  var lastE = this.Field(x,y);

  var traceFinished = false;
  var nstep = 0;
  while(true) {
    nstep++;
    var E = this.Field(x,y);
    
    var dx = E.gx * step *fieldline.dir;
    var dy = E.gy * step *fieldline.dir;

    // Parasitic calculation. Find line segments that cross equipotential lines.
    if(this.do_equipotential) 
    {
      var span = SpansIntegerMultiple(lastE.U, E.U, potential_multiple);
      if(span!=null) {
        pnode = { U: span*potential_multiple, E1: lastE, E2: E };
        this.potentialnodes.push(pnode);
      }
    }
    
    // This check doesn't work.
    // I've also tried the form of E/U, which should check for small field in high-potential areas (which is what I want to check)
    // I can't just check for small field, because that happens naturally at large distances from the middle.
    // if(Math.abs(E.E) < 1 && Math.abs(E.U) > 2) {
    //   console.log("Line went through effective zero-field region. This isn't a good line. E=",E.E," U=",E.U, " step=",step);
    //   return false;
    // }
    
    x += dx;
    y += dy;
    fieldline.points.push({x:x,y:y});        
    lastE = E;

  
    var collide = this.FindCollision(x,y);
    if(collide && (fieldline.dir*collide.q < 0) && nstep>1) {
      // Find the best possible node for this line.
      if(collide.n_nodes >= collide.nodes.length+1 == 0) {
        // Comment these lines out if you want it to just sail through without stopping...
        // console.warn("Line failed - hit q=",collide.q,"which has no nodes left.");
        // return false; //nodeFinished=false;
      } else {
        this.DoCollision(collide,x,y);
        fieldline.endCharge = collide;
        fieldline.nstep = nstep;
        console.log("Line succeeded - hit q=",collide.q);
        return true; // nodeFinished
      }
    }
              
    if(nstep>max_steps){
      fieldline.endCharge = null;
      fieldline.endAngle     = null;
      fieldline.endNodeAngle = null;
      fieldline.nstep = nstep;
      console.log("Line succeeded - no hit");
      return true;
    }  // if nstep 
  } // trace loop
}
Applet.prototype.FindFieldLines = function(){
  

  this.fieldLines = [];
  this.potentialnodes = []; 
  this.equipotential_lines = [];


  var total_charge = 0;
  var max_x = -1e20;
  var min_x =  1e20;
  var max_y = -1e20;
  var min_y =  1e20;
  for(var i=0 ;i<this.charges.length; i++) {
    var charge = this.charges[i];
    total_charge += charge.q;
    charge.r = 0.12*Math.sqrt(Math.abs(charge.q));
    charge.n_nodes = Math.round(Math.abs(source_lines_per_unit_charge*charge.q));
    charge.nodes = []; // All successful or unsuccesful nodes
    charge.nodesUsed = []; // Nodes that have actually worked.
    charge.nodesNeeded = []; // Some idea what nodes we should try.
    if(charge.x > max_x) max_x = charge.x;
    if(charge.x < min_x) min_x = charge.x;
    if(charge.y > max_y) max_y = charge.y;
    if(charge.y < min_y) min_y = charge.y;
  }
  

  // rank them. Use minority charge carriers first: their nodes HAVE to connect.
  this.charges.sort(chargesort);
  if(total_charge<0) this.charges.reverse();

  console.log("Doing escaping lines -------------- ");
  // Find fieldlines that come from outside the area, assuming there is a majority charge carrier.
  var escaping_lines = Math.abs(total_charge* source_lines_per_unit_charge);
  for(var i=0;i<escaping_lines;i++) {
    console.log("Doing escaping line.");
    // Find a position very far away from the charges.
    var r = Math.max(this.xmax,this.ymax) * 10;
    if(isNaN(r)) r = 10;
    var theta = i*TWO_PI/escaping_lines;
    var x =  r*Math.cos(theta);
    var y =  r*Math.sin(theta);
    
    var fieldline = { startCharge: null };
    if(total_charge > 0)  fieldline.dir = -1;
    else                  fieldline.dir = 1;
    fieldline.start_x = x;
    fieldline.start_y = y;
    fieldline.start = "outside";
    var nodeFinished = this.TraceFieldLine(fieldline); 
    if(nodeFinished) {
      this.fieldLines.push(fieldline);      
    } else {
      console.log("incoming line failed");
    }
    
  }
  


  // Now loop through again, finding unused nodes and tracing field lines from those
  // nodes until they either hit another charge or they require too many computational cycles.

  for(var i=0 ;i<this.charges.length; i++) {
    var random_seed = 0;
    var charge = this.charges[i];    
    // console.log("Find field lines for charge ",i," with charge ",charge.q);
    this.ctx.fillStyle = 'blue';
    console.log("Doing charge",i,"with q=",charge.q,"which has ",charge.nodesUsed.length,"/",charge.n_nodes," nodes");


    while(charge.nodesUsed.length < charge.n_nodes && charge.nodes.length<source_lines_per_unit_charge*5) {
      if(charge.nodes.length>source_lines_per_unit_charge*4) {
        console.warn("Wow! Tried way too many nodes.",charge.nodes);
      }
      console.log("Doing node on charge",i);
      
      var start_angle = this.FindNodePosition(charge);

      var r = charge.r;
      // Boost in initial direction by radius.
      var fieldline = { startCharge: charge };
      fieldline.start = "charge";
      var nodeFinished = false;

      // console.log("Try: ",nodeTries,"Trying angle:",start_angle*180/Math.PI,nodeTries);
      fieldline.start_x = charge.x + charge.r* Math.cos(start_angle);
      fieldline.start_y = charge.y + charge.r* Math.sin(start_angle);
      fieldline.start_angle = start_angle;
      var dir = 1;
      if(charge.q<0) dir = -1;
      fieldline.dir     = dir;
      
      var nodeFinished = this.TraceFieldLine(fieldline); 
      if(nodeFinished) {
        charge.nodesUsed.push(start_angle);
        this.fieldLines.push(fieldline);      
      }
    } // nodeFinished
  }
  
  if(this.do_equipotential){
  // Find equipotential lines.
  // Trace around all the equpotenial nodes we've found.
  console.log("looking at potentialnodes: ", this.potentialnodes.length);

  while(this.potentialnodes.length>0) {
    var pnode = this.potentialnodes.shift();
    var Utarget = pnode.U;
    // Fresh node. Approximate the point of best potential.
    console.log("Trying node, Utarget=",Utarget);

    var point = pnode.E1;  // Pick one of the two end segments.
    point = this.FindPositionOfU(point,Utarget,Utolerance);
    var xstart = point.x;
    var ystart = point.y;
    // var fU = (pnode.U - pnode.E1.U)/(pnode.E2.U-pnode.E1.U);
    // if(pnode.E2.U < pnode.E1.U) fU = -fU;
    // var startx = pnode.x1 + fU*(pnode.x2-pnode.x1);
    // var starty = pnode.y1 + fU*(pnode.y2-pnode.y1);
    // var startE = this.Field(startx,starty);
    // Start tracing this equpotential back and forth.
    for(var dir = -1; dir<3; dir +=2) {
      var line = {U: Utarget, points:[{x:point.x,y:point.y}]};
      var done = false;
      // console.log("start line at",startE.U,pnode.U,pnode);
      var np = 0;
      while(!done) {
        np++;
        // console.log(point);
        var newx =  point.x + point.gy * step_equi * dir;  // Not a typo. .
        var newy =  point.y - point.gx * step_equi * dir;  // We're going perpendicular to the field!
        var next_point = this.Field(newx,newy);        
        var next_point = this.FindPositionOfU(next_point,Utarget,Utolerance); // refine
        
        // Check for intersection with other potentialnodes. Delete them as we go.
        for(var i=0;i<this.potentialnodes.length;i++) {
            var othernode = this.potentialnodes[i];
            if(othernode.U == Utarget) {
              if(LineSegmentsIntersect(point,next_point,
                                       othernode.E1, othernode.E2)) {
              console.warn("collide with node!  left:",this.potentialnodes.length);
              this.potentialnodes.splice(i,1);
            }
          }
        }
        // var d2 = (next_point.x - xstart)*(next_point.x - xstart)+(next_point.y - ystart)*(next_point.y - ystart);
        // console.log("distance from start: ",Math.sqrt(d2));
        if(np>2 && LineSegmentsIntersect(point,next_point,
                                          pnode.E1,pnode.E2))  {
          done = true;
          dir+=2; // exit dir loop
          console.warn("looped line");
        } else if(np>max_equi_step){
          console.warn('gave up on line');
          done = true;
        } 
        line.points.push({x:next_point.x,y:next_point.y});
        point = next_point;
        // console.log(E.U);
      }
      this.equipotential_lines.push(line);
      console.log("End U at",point.E.U);
    }
    // break;
  }
  }
 
  
}
Applet.prototype.TotalEnergy = function(){
  var tot = 0;
  for(var i=1 ;i<this.charges.length; i++) {
    for(var j=0 ;j<i; j++) {
    
      // compute potential at this point for all other charges.
      var ci = this.charges[i];
      var cj = this.charges[j];

      var dx = ci.x-cj.x;
      var dy = ci.y-cj.y;
      var r2 = dx*dx+dy*dy;
      var r = Math.sqrt(r2);
      tot += ci.q*cj.q/r;
    }
  }
  return tot;
}
Applet.prototype.Draw = function(){
  this.Clear();
  this.ctx.save();
  
  this.do_equipotential = false;
  // uncomment this to is you want it to make equipotentials, only worked when it was in 2d
  // this.do_equipotential = $('#ctl-do-eqipotential').is(":checked");
  
  this.canvas_translate = { x: this.canvas.width/2, y: this.canvas.height/2};
  this.canvas_scale     = { x: this.canvas.width/this.width_x, y: -this.canvas.width/this.width_x};
  
  this.ctx.translate(this.canvas_translate.x,this.canvas_translate.y);
  this.ctx.scale(this.canvas_scale.x,this.canvas_scale.y);
  this.xmin = -this.width_x/2;
  this.xmax =  this.width_x/2;
  this.ymin = -this.width_x/2 * this.canvas.height/this.canvas.width;
  this.ymax =  this.width_x/2 * this.canvas.height/this.canvas.width;
  
  this.ctx.strokeStyle = 'black';
  this.ctx.lineWidth = 0.01;
  this.ctx.beginPath();
  this.ctx.moveTo(this.xmin,this.ymin);
  this.ctx.lineTo(this.xmax,this.ymin);
  this.ctx.lineTo(this.xmax,this.ymax);
  this.ctx.lineTo(this.xmin,this.ymax);
  this.ctx.lineTo(this.xmin,this.ymin);
  this.ctx.stroke();
  this.ctx.beginPath();
  
  var urlparams = "";
  for(var i = 0; i< this.charges.length; i++) {
    if(i==0) urlparams+="?";
    else   urlparams += "&";
    urlparams += "q"+i+"=";
    urlparams += this.charges[i].q + "," + parseFloat(this.charges[i].x.toFixed(3)) + "," + (parseFloat(this.charges[i].y.toFixed(3)));
    if(hide_charge_values) urlparams += "&hideq=1";
    urlparams+="&lines=" + source_lines_per_unit_charge;
  }
  $('#totalenergy').html("Total Energy: "+this.TotalEnergy().toFixed(1));
  $('#linktothis').attr('href',urlparams);

  console.time("FindFieldLines");
  this.FindFieldLines()
  console.timeEnd("FindFieldLines");
  
  this.DrawFieldLines();
  this.DrawCharges();
  if(this.do_equipotential) this.DrawEquipotentialLines();
  
  this.ctx.restore();
  
}
Applet.prototype.DrawFieldLines = function(){
  console.time("Drawing lines")
  this.ctx.lineWidth = 0.02;
  for(var i=0;i<this.fieldLines.length;i++) {
    var line = this.fieldLines[i];
    //console.log("Drawing line ",i);
    this.ctx.strokeStyle = 'black';
    // this.ctx.strokeStyle = 'blue';
    // if(line.startCharge.q >0) this.ctx.strokeStyle = 'red';
    this.ctx.beginPath();
    this.ctx.lineJoin = "round";
    this.ctx.moveTo(line.start_x,line.start_y)
    for(var j=1;j<line.points.length;j++) {
      var p = line.points[j];
      this.ctx.lineTo(p.x,p.y);
    }
    this.ctx.stroke();
    
	
	
    var n = line.points.length;
    // Add arrow. Find the midway point along the line.
    var j = Math.round((n-1)/2);
    // console.log(j,line.points.length);
    var x = line.points[j].x;
    var y = line.points[j].y;
    // Ensure arrow is on the screen - keep halving the midway point until we reach it.
    while(x<this.xmin || x>this.xmax || y<this.ymin || y>this.ymax) {
      if(line.start == "outside") j = Math.round(n-(n-j)/2);
      else                        j = Math.round(j/2);
      x = line.points[j].x;
      y = line.points[j].y;
      //console.log(j);
      if(j<=1 || j>=n-3) break;
    }
    dx = line.dir*(line.points[j+1].x - x);
    dy = line.dir*(line.points[j+1].y - y);
    this.ctx.save();
    this.ctx.translate(x,y);
    this.ctx.fillStyle = 'black';
    var angle = (Math.atan2(dy,dx)+TWO_PI)%TWO_PI;
    this.ctx.rotate(angle);
    var lx = 0.2;
    var ly = 0.08;
    this.ctx.beginPath();
    this.ctx.moveTo(lx,0);
    this.ctx.lineTo(0,ly);
    this.ctx.lineTo(0,-ly);
    this.ctx.lineTo(lx,0);
    this.ctx.closePath();
    this.ctx.fill();
    this.ctx.restore();


    // Make dots.
    // for(var j=0;j<line.points.length;j++) {
    //   this.ctx.beginPath();
    //   var p = line.points[j];
    //   this.ctx.arc(p.x,p.y,0.01,0,Math.PI*2,true);
    //   this.ctx.fill(); 
    // }
  }
  console.timeEnd("Drawing lines")

}
Applet.prototype.DrawEquipotentialLines = function(){
  console.time("Drawing potential lines")

  for(var i=0;i<this.equipotential_lines.length;i++) {
    var line = this.equipotential_lines[i];
    this.ctx.beginPath();
    this.ctx.lineWidth = 0.02;
    this.ctx.strokeStyle = "green";
    this.ctx.lineJoin = "round";
    this.ctx.moveTo(line.points[0].x,line.points[0].y)
    for(var j=1;j<line.points.length;j++) {
      var p = line.points[j];
      this.ctx.lineTo(p.x,p.y);
    }
    this.ctx.stroke();

    this.ctx.strokeStyle = "black";
    this.ctx.lineWidth = 0.01;
    
    // Make dots.
    // for(var j=0;j<line.points.length;j++) {
    //   this.ctx.beginPath();
    //   var p = line.points[j];
    //   this.ctx.arc(p.x,p.y,0.01,0,Math.PI*2,true);
    //   this.ctx.stroke();
    // }
  }
  console.timeEnd("Drawing potential lines")

}
Applet.prototype.DrawCharges = function(){
  // Draw charges. Do this last so line tails are covered.
  for(var i=0 ;i<this.charges.length; i++) {
    var charge = this.charges[i];    
    this.ctx.fillStyle = 'blue';
    if(charge.q >0) this.ctx.fillStyle = 'red';
    if(charge.highlight) this.ctx.lineWidth = 0.03;
    else                 this.ctx.lineWidth = 0.01;
    var x = charge.x;
    var y = charge.y;
    var r = charge.r;
    //console.log(charge.x,charge.y,charge.r,0,Math.PI*2,true);
    this.ctx.beginPath();
    this.ctx.arc(x,y,r,0,Math.PI*2,true);
    this.ctx.fill();
    this.ctx.stroke();

    // Draw attempted node positions, for debugging
    // for(var j=0;j<charge.nodes.length;j++) {
    //   var t= charge.nodes[j];
    //   x = charge.x + r*Math.cos(t);
    //   y = charge.y + r*Math.sin(t);
    //   this.ctx.beginPath();
    //   this.ctx.arc(x,y,r/5,0,Math.PI*2,true);
    //   this.ctx.stroke();
    // }

    this.ctx.save();
    this.ctx.translate(charge.x,charge.y);
    this.ctx.scale(0.01,-0.01);
    this.ctx.fillStyle = 'white';
    this.ctx.strokeStyle = 'white';
    this.ctx.textBaseline = 'middle';
    this.ctx.textAlign = 'center';
    this.ctx.font = '12pt sans-serif';
    var s;
    if(charge.q<0) s = "-";
    else           s = "+";
    s += parseInt(Math.abs(charge.q));
    if(!hide_charge_values) // protect against old browsers
      this.ctx.fillText(s,0,0);
    this.ctx.restore();
  }
  
}
function getAbsolutePosition(element) {
   var r = { x: element.offsetLeft, y: element.offsetTop };
   if (element.offsetParent) {
     var tmp = getAbsolutePosition(element.offsetParent);
     r.x += tmp.x;
     r.y += tmp.y;
   }
   return r;
 };
 Applet.prototype.GetEventXY = function(ev){
  // Convert mouse click coordinates to the mathematical plane.
   var offset = getAbsolutePosition(this.canvas);
   var x = ev.pageX;
   var y = ev.pageY;
   //$('#debug').html("DoMouse "+ ev.type + " " + ev.originalEvent.touches.length + " " + x +  " " + y);    
  
   if((ev.type =='touchstart') || (ev.type =='touchmove') || (ev.type =='touchend')) {
     ev.preventDefault();
     //$('#debug').html("DoMouse "+ ev.type + " " + ev.originalEvent.touches.length + " " + x +  " " + y);    
     x = ev.originalEvent.touches[0].pageX;
     y = ev.originalEvent.touches[0].pageY;
   }
   x = x - offset.x;
   y = y - offset.y;    
   x -= this.canvas_translate.x;
   y -= this.canvas_translate.y;
   x /= this.canvas_scale.x;
   y /= this.canvas_scale.y;
   return {x:x, y:y};
 }
 Applet.prototype.DoMouse = function(ev){
  var xy = this.GetEventXY(ev);
  var x = xy.x;
  var y = xy.y;

  // console.log(ev.type,x,y);

  var update = false;
  

  if(ev.type === 'mousedown' || ev.type ==='touchstart') {
    // See if we're clicking a charge.
    var charge = this.FindCollision(x,y);
    if(charge) {
      this.dragging = true;
      this.charge_dragged = charge;
      charge.highlight = true;
      update = true;
    }
  }
  if(ev.type === 'mousemove' || ev.type ==='touchmove') {
    if(this.dragging) {
      this.charge_dragged.x = x;
      this.charge_dragged.y = y;
      update = true;    
    }
  }
  if(ev.type === 'mouseup' || ev.type ==='touchend') {
    if(this.charge_dragged) this.charge_dragged.highlight = false;
    this.charge_dragged = null
    this.dragging = false;
    update = true;
  }
  
  if(ev.type === 'mouseout') {
    if(this.charge_dragged){
      // find it in the list.
      var which = 0;
      for(var i=0;i<this.charges.length;i++) if(this.charge_dragged==this.charges[i]) which=i;
      this.charges.splice(which,1);
      this.charge_dragged = false;
      this.dragging = false;
      update = true;
    }
  }
    
  if(update) this.Draw();
}
Applet.prototype.AddCharge = function(ev){
  console.log("AddCharge",ev);
  var q = parseFloat(ev.currentTarget.getAttribute('q'));
  var xy = this.GetEventXY(ev);
  var x = xy.x;
  var y = xy.y;
  
  var charge = { q : q,  x : x,  y: y , r: 0.12*Math.sqrt(Math.abs(q))};
  this.charges.push(charge);
  
  this.dragging = true;
  this.charge_dragged = charge;
  charge.highlight = true;
  update = true;
  
  this.Draw();
}
Applet.prototype.AddChargeRandom = function(ev){
  console.log(ev);
  var q = parseFloat(ev.currentTarget.getAttribute('q'));
  console.log(q);
  this.xmin = -this.width_x/2;
  this.xmax =  this.width_x/2;
  this.ymin = -this.width_x/2 * this.canvas.height/this.canvas.width;
  this.ymax =  this.width_x/2 * this.canvas.height/this.canvas.width;
  var x = (Math.random()*1.8 - 0.9)*(this.xmax-this.xmin)/2;
  var y = (Math.random()*1.8-0.9) *(this.ymax-this.ymin)/2;
  this.charges.push({
    q : q,  x : x,  y: y , r:0.12*Math.abs(q)
  });

  this.Draw();
}
Applet.prototype.PrintHQ = function(){
  console.log("Applet::PrintHQ");
  // First, save our current state.
  var saveCanvas = this.canvas;
  var saveCtx    = this.ctx;
  var saveWidth  = this.width;
  var saveHeight = this.height;

  // Second, create an offscreen drawing context.
  var canvas = document.createElement("canvas");
  this.canvas = canvas;
  canvas.width = saveWidth * gPrintScale;
  canvas.height = saveHeight * gPrintScale;
  // this.width  = saveWidth * gPrintScale;
  // this.height = saveHeight * gPrintScale;
  this.ctx = this.canvas.getContext("2d");

  // Now do the actual draw
  // this.Resize();
  this.ctx.scale(gPrintScale,gPrintScale);
  this.Draw();

  // Save the print result
  gPrintBuffer = this.canvas.toDataURL('image/png');


  // Restore defaults.
  this.canvas = saveCanvas;
  this.ctx    = saveCtx;


  // nuke buffer.
  // document.removeChild(canvas); // Doesn't work. Is this thing a memory leak? I don't think so - I htink the canvas vanishes when it goes out of scope.
};
*/