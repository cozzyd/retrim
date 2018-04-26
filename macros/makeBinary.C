void makeBinary(const char * file) 
{

  gSystem->Load("lib/libretrim.so"); 

  retrim::CollisionReader c(file); 
  TString newname =  TString(file).ReplaceAll(".txt",".coll"); 
  cout << "writing " << newname << endl; 
  c.writeBinary(newname); 
}
