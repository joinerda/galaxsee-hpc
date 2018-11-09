#macro basicView(lookFrom,lookTo,Angle)
    #local axisRadius = 0.05;
    background { color rgb<0.9,0.9,0.9> }
    camera {
      location lookFrom
      look_at lookTo
      right<-1.33,0,0>
      angle Angle
    }      
    sky_sphere {pigment {color rgb<0,0,0>}}

    global_settings { ambient_light rgb<5, 5, 5> }   
//    light_source { <2, 4, -3> color White}
    #local light = lookFrom + (lookFrom-lookTo)*0.1;
    light_source { light color White}
/*
    cylinder {
        <0,0,0>,<4,0,0>,axisRadius
        texture {
            pigment { color Red }
        }
    }
    cylinder {
        <0,0,0>,<0,4,0>,axisRadius
        texture {
            pigment { color Green }
        }
    }
    cylinder {
        <0,0,0>,<0,0,4>,axisRadius
        texture {
            pigment { color Blue }
        }
    }
*/
    #local i=-10;
    #local gridBright = 0.4;
    #local gridTrans = 0.0;
/*
    #while (i<10)
        cylinder {
            <i,0,-10>,<i,0,10>,0.25*axisRadius
            texture {
                pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
            }
        }
        cylinder {
            <-10,0,i>,<10,0,i>,0.25*axisRadius
            texture {
                pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
            }
        }
        #local i = i + 1;
    #end
*/
    cylinder {
        <-1,-1,-1>,<-1,-1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,-1,-1>,<1,-1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,1,-1>,<1,1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <-1,1,-1>,<-1,1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }

    cylinder {
        <1,-1,-1>,<-1,-1,-1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,-1,1>,<-1,-1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,1,-1>,<-1,1,-1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,1,1>,<-1,1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }

    cylinder {
        <-1,-1,-1>,<-1,1,-1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,-1,-1>,<1,1,-1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <-1,-1,1>,<-1,1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
    cylinder {
        <1,-1,1>,<1,1,1>,0.25*axisRadius
        texture {
            pigment { color rgbt<gridBright,gridBright,gridBright,gridTrans> }
        }
    }
#end

