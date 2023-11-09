// This file is a wrapper for a call to MathJax TeX/Latex formula 
// rendering system. If you have MathJax installed locally, you may
// want to modify this file so that "filename" points to the local 
// version of MathJax.js .  In this case you will not need an 
// internet connection while reading "Medel Description"


  // set remote location of MathJax.js script
  var filename="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

  // or, if it is installed locally, set a local path to it
  //  var filename="/opt/mathjax/mathjax-MathJax-07669ac/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

  var head = document.getElementsByTagName('head')[0];
  var script = document.createElement('script');
  script.type = 'text/javascript';
  script.onload = function(){
    // remote script has loaded
    script.onload = null;
    MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
    };

  script.src = filename;

  head.appendChild(script);
