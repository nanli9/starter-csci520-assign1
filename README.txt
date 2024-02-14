<Please submit this file with your solution.>

CSCI 520, Assignment 1

<Nan Li>

================
Please run the program in Windows machine, the stencil selection is not working on my mac
<Description of what you have accomplished>
Description:
Finish the homework following the homework instructions. First, the structural, shear, and bend springs informations are stored for each point.
Then, external force field is applied to the jello cube by trilinear interpolation, since the force field is given in grid format. After that, penalty
method is applied to implement the collision. The collision detection is determined by whether f(x,y,z) < 0, all the normal vectors for each plane are
set inward to the center of the origin.


<Also, explain any extra credit that you have implemented.>
Extra credit:
1.add an incline plane with collision detection and handler. Since there is no incline plane data in given world files, I manually add it in jello.cpp line 339. if you want to read
from a world file, just comment it out.
2.add an new interactive key "t": it will load a texture into the bounding box, the texture will be discarded if "t" is pressed again.
(stb_image.h file is added to load the image."https://github.com/nothings/stb/blob/master/stb_image.h")
3.Users can apply force to the entire jello by dragging the mouse at any place of the screen.
4. use key 'q' to toggle whether users can select the single point then act the force on it(wireframe mode only). After press 'q', users can
move the mouse to the point and the point will become red if selected. Then users can apply force to the single point by using the left mouse key
to drag. You can change coefficient for user input force in input.cpp line 68 and line 83.