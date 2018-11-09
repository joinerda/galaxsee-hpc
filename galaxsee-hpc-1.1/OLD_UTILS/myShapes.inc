#include "colors.inc"    
#include "myTextures.inc"

#macro arrow(Base,Tip,aScale,Color)
    #local Trans=0.1;
    #local Break = Base+(Tip-Base)*0.6;
    #local Distance = vdot(Base-Break,Base-Break);
    sphere {
      Base, min(Distance,0.4)*aScale
      texture {
        pigment { color rgbt<Color.red,Color.green,Color.blue,Trans> }
      }
    }
    #if(Distance > 0.01)
    #local CylColor = (Color+White)*0.5;
        cylinder {
          Base,Break,min(Distance,0.3)*aScale
          texture {
            pigment { color
                rgbt<CylColor.red,CylColor.green,CylColor.blue,Trans> }
          }
        }
    #end
    #local Distance = vdot(Base-Break,Base-Break);
    #if(Distance > 0.01)
        cone {
          Break,min(Distance,0.6)*aScale,Tip,0.0
          texture {
            pigment { color rgbt<1,1,1,Trans> }
          }
        }
    #end
#end

#macro point(position,pScale,Trans,Color)
    sphere {
      position, 0.4*pScale
      texture {
        pigment { color rgbt<Color.red,Color.green,Color.blue,Trans> }
      }
    }
#end

#macro metal_bar(Position1,Position2,Radius)
cylinder {
    Position1,Position2,Radius
    texture {
        pigment { color
            rgbt<0.6,0.6,0.9,0.15>
        }
    }
}
#end


