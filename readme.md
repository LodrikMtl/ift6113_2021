# IFT 6113 (Fall 2021)
Template repository for homeworks.


### Usage
If you are new to git (or github), the easiest way to start up with this code is to just download it as zip-folder. Find this option on tope (green button "code"). Alternatively, you can fork this repo, and then update it as new homework templates arrive. 

## HW instructions
See instructions in subfolders.

***

## Homework 1
Your task is to develop subdivision schemes. 14 Sept - 1 Oct

### Instructions to run the given codes
If you are on Windows, you can copy the example.* project two times and rename one butterflyloop and the other one sqrt3. Add mainlb.cpp to butterflyloop and mains.cpp to sqrt3. Set the project needed to be run as startup project and launch it.

I struggle to change the cmake file, so if it is not fine, I can try to submit a version that will setup correctly the code.

#### Loop+Butterfly Arguments
mainlb.cpp expect the arguments in the following orders: ```{method_name} {n} {input_path}``` where method_name $$\in ["butterfly","loop"]$$.

#### Sqrt(3) Arguments
mains.cpp expect the arguments in the following orders: ```{n} {input_path} ```.
