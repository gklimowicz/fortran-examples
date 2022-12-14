import sys

class RundeckSectionOld:
    def __init__(self,
                 name="",
                 components=[],
                 obj_files=[],
                 input_files={},
                 parameters_text="XX",
                 parameters={},
                 variants={}):
        #print "name>>> " + name
        #print "obj_files>>> ", obj_files
        # use : self.obj_files = list(obj_files)
        #       self.input_files = dict(input_files)
        self.name = name
        self.obj_files = []
        self.obj_files = obj_files
        self.components = components
        self.input_files = input_files
        self.parameters_text = parameters_text
        #print ">>> ", self.parameters_text
        self.parameters = parameters
        self.opt = []
        self.variants = variants


    def update(self, x):
        self.name += x.name
        self.obj_files.extend(x.obj_files)
        self.components.extend(x.components)
        self.input_files.update(x.input_files)
        self.parameters_text += x.parameters_text
        self.parameters.update(x.parameters)
        self.opt.extend(x.opt)
        self.variants.update(x.variants)


    def prt(self):
        print
        print "##################################################"
        print "#Section: "+self.name
        print
        print "Object modules:"
        for name in self.obj_files:
            print name,

        print "\n"
        print "Components:"
        for name in self.components:
            print name,

        print
        print "Component Options:"
        for i in self.opt:
            print i,

        print "\n"
        print "Data input files:"
        for i in self.input_files.keys():
            print i+"="+self.input_files[i]

        print
        print "&&PARAMETERS"
        print self.parameters_text
        for i in  self.parameters.keys():
            print i+"=",self.parameters[i]
            if isinstance(self.parameters[i],list):
                print str(self.parameters[i])[1:-1]

    
    def __call__(self,*varlist):
        x = RundeckSection()
        x.update(self)
        #x.name = self.name
        #x.obj_files.extend(self.obj_files)
        #x.components.extend(self.components)
        for var in varlist:
            x.update(self.variants[var])
            #x.obj_files.extend(self.variants[var].obj_files)
            #x.parameters.update(self.variants[var].parameters)
        return x



def f(**kwargs):
    x = dict()
    for k, v in kwargs.iteritems():
        x[k] = v
    return x

# functions for RS constructor

def name(x):
    return("name",x)

def components(*args):
    return("components", args)

def component_options(*args):
    return("component_options", args)

def obj_files_text(x):
    return("obj_files_text",x)

def obj_files(*args):
    return("obj_files", args)

def input_files_text(x):
    return("input_files_text",x)

def input_files(**kwargs):
    x = dict()
    for k, v in kwargs.iteritems():
        x[k] = v
    return("input_files", x)

def parameters_text(x):
    return("parameters_text",x)

def parameters(**kwargs):
    x = dict()
    for k, v in kwargs.iteritems():
        x[k] = v
    return("parameters", x)

def variants(**kwargs):
    x = dict()
    for k, v in kwargs.iteritems():
        x[k] = v
    return("variants", x)



class RundeckSection:
    def __init__(self,*args):
        self.name = ""
        self.components = list()
        self.component_options =list()
        self.obj_files_text = ""
        self.obj_files = list()
        self.input_files_text = ""
        self.input_files = dict()
        self.parameters_text = ""
        self.parameters = dict()
        self.variants = dict()
        # kwargs is a dictionary.
        for i in args:
            #print i[0],"--->",i[1]
            if isinstance(i[1],dict):
                self.__dict__[i[0]].update(i[1])
            elif isinstance(i[1],list):
                self.__dict__[i[0]].extend(i[1])
            else:
                self.__dict__[i[0]] += i[1]

    def update(self, x):
        self.name += x.name
        self.components.extend(x.components)
        self.component_options.extend(x.component_options)
        self.obj_files_text += x.obj_files_text
        self.obj_files.extend(x.obj_files)
        self.input_files_text += x.input_files_text
        self.input_files.update(x.input_files)
        self.parameters_text += x.parameters_text
        self.parameters.update(x.parameters)
        self.variants.update(x.variants)

    def __call__(self,*varlist):
        x = RundeckSection()
        x.update(self)
        for var in varlist:
            x.update(self.variants[var])
        return x

    def prt(self):
        print
        print "##################################################"
        print "#Section: "+self.name
        print
        print "Object modules:"
        for name in self.obj_files:
            print name,

        print "\n"
        print "Components:"
        for name in self.components:
            print name,

        print
        print "Component Options:"
        for i in self.component_options:
            print i,

        print "\n"
        print "Data input files:"
        for i in self.input_files.keys():
            print i+"="+self.input_files[i]

        print
        print "&&PARAMETERS"
        if self.parameters_text != "": print self.parameters_text
        for i in  self.parameters.keys():
            print i+"=",self.parameters[i]
            if isinstance(self.parameters[i],list):
                print str(self.parameters[i])[1:-1]


    def prt_obj_files(self):
        print "!------------",self.name,"------------"
        print self.obj_files_text,
        for name in self.obj_files:
            print name,
        print

    def prt_components(self):
        if not self.components : return
        print "!------------",self.name,"------------"
        for name in self.components:
            print name,
        print

    def prt_component_options(self):
        if not self.component_options : return
        print "!------------",self.name,"------------"
        for name in self.component_options:
            print name,
        print

    def prt_input_files(self,res):
        if not self.input_files and not self.input_files_text: return
        if res == "M20":
            print "! M20"
            r = {"$resij":"72x46","$resIJ":"72X46",
                 "$resdeg":"4x5","$resDEG":"4X5"}
        elif res == "F40":
            print "! F40"
            r = {"$resij":"144x90","$resIJ":"144X90",
                 "$resdeg":"2x2.5","$resDEG":"2X2.5"}
        else:
            sys.exit("Unknown resolution")
        
        print "!------------",self.name,"------------"
        s = self.input_files_text
        #print "processing", s 
        for k in r.keys():
            #print "trying:",k,r[k]
            s = s.replace(k,r[k])
        print s
        for i in self.input_files.keys():
            s = self.input_files[i]
            #print "processing", s 
            for k in r.keys():
                #print "trying:",k,r[k]
                s = s.replace(k,r[k])
            print i+"="+s

    def prt_parameters(self):
        if not self.parameters and not self.parameters_text: return
        print "!------------",self.name,"------------"
        print self.parameters_text
        for i in self.parameters.keys():
            print i+"="+self.parameters[i]
            if isinstance(self.parameters[i],list):
                print str(self.parameters[i])[1:-1]



def prt_rundeck(title,descr,cppopts,sect,namelist,res):
    print title
    print
    print descr
    print
    print "Preprocessor Options"
    for i in cppopts:
        print i
    print "End Preprocessor Options"
    print

    # Sections:
    print "Object modules:"
    for i in sect:
        i.prt_obj_files()
    print

    print "Components:"
    for i in sect:
        i.prt_components()
    print

    print "Component Options:"
    for i in sect:
        i.prt_component_options()
    print

    print "Data input files:"
    for i in sect:
        i.prt_input_files(res)
    print

    print "Label and Namelist:  (next 2 lines)"
    print title
    print

    print "&&PARAMETERS"
    for i in sect:
        i.prt_parameters()
    print "&&END_PARAMETERS"

    print namelist

