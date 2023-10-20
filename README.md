# Satellite-Navigation

**The application allows you to determine DOP (Dilutions of precision) coefficients and satellite positions for any selected day and place on earth using almanac file.**

In file main.py you can change the inputs - start date, end date, how often you need to calculate the results, position of the receiver, mask (how much above the horizon must be satelite to take it into acount during calculation), and numbers of the satellites you want to take into considaration (if you left this table empty, the porgram will calculate everything for any satelite visible from given place).

## Output:

**Plot of the values of the DOP coefficients in time (in hours) and the visibility of the satellites in time (taking into consideration given mask).**

![image](https://github.com/MariaMank/Satellite-Navigation/assets/92314221/48bd088a-1992-4d23-8f60-48d72b1957ff)

**The diagram of the elevation of the given satelites.**

![image](https://github.com/MariaMank/Satellite-Navigation/assets/92314221/fc1f4f07-85e1-4053-a579-d1c1e47ebd07)

**Aimation window which presents the visibility and position of satellites**

![image](https://github.com/MariaMank/Satellite-Navigation/assets/92314221/7ab7f8f8-7890-4ac5-a1a3-7bc9a0bbd4fa)

I also uploaded the video presenting how does the animation work - 2022-04-08 01-08-50.mp4

File main.py is fully written by me and contains all of the most important in this project computations.
Files date2tow.py, groundstack_stud.py, read_yuma.py was created by Maciej Grzymała - teacher of this course.
File skyplot2.py was also created by Maciej, but I updated it to make animation.

