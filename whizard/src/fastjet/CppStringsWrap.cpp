#include<string>

using namespace std;

extern "C" {

  void cpp_str_delete (string* str) {
    delete str;
  }

  int cpp_str_length (const string* str) {
    return str->length ();
  }

  char cpp_str_get (const string* str, const int i) {
    return (*str)[i];
  }

}
