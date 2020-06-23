#ifndef PROJECT_MY_EXCEPTION_H
#define PROJECT_MY_EXCEPTION_H

#include "iostream"
#include "exception"

using namespace std;

class my_exception : public exception {
public:
    explicit my_exception(const char *errMsg = "a") : _errMsg(errMsg){}

public:
    const char *what() const throw() {

        cout<<"Exception caught!\nError msg:"<<_errMsg.c_str()<<endl;
        return "C++ Exception\n";
    }

private:
    string _errMsg;
};

#endif //PROJECT_MY_EXCEPTION_H
