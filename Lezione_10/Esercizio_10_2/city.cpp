#include "city.h"

using namespace std;

void city :: initialize(){
   _pos.resize(_ndim);
   _pos_old.resize(_ndim);
   return;
}

vec city :: get_pos(){
    return _pos;
}

void city :: set_pos(vec pos){
    _pos = pos;
}

vec city :: get_pos_old(){
    return _pos_old;
}

void city :: set_pos_old(vec pos_old){
    _pos_old = pos_old;
}

string city :: get_name(){
    return name;
}

void city :: set_name(string name_new){
    name = name_new;
}
