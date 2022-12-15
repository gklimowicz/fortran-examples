-- Copyright (c) 2012 Harald Klimach <harald@klimachs.de>
--
-- Parts of this file were written by Harald Klimach for German Resarch
-- School of Simulation Sciences.
--
-- Permission is hereby granted, free of charge, to any person obtaining a copy
-- of this software and associated documentation files (the "Software"), to deal
-- in the Software without restriction, including without limitation the rights
-- to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
-- copies of the Software, and to permit persons to whom the Software is
-- furnished to do so, subject to the following conditions:
--
-- The above copyright notice and this permission notice shall be included in
-- all copies or substantial portions of the Software.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
-- IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
-- FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
-- IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
-- DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
-- OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
-- OR OTHER DEALINGS IN THE SOFTWARE.
-- ************************************************************************** --

-- Small sample Config File in Lua
width = 200
height = 150

stl_files = { {'filename', 123, 'binary'},
              {'geomfile', 456, fileformat='ascii'} }

coord = { 0.1, 0.0, 0.0 }

--function 'gausspulse'
function gauss_pulse(x, y, z, origin, amplitude, hwidth)
  fact = -0.5/(hwidth*hwidth)
  dist = (x - origin[1])*(x - origin[1]) 
       + (y - origin[2])*(y - origin[2]) 
       + (z - origin[3])*(z - origin[3])
  res = amplitude * math.exp(fact*dist)
  return res
end

--function 'ic_density' 
function ic_density(x, y, z)
  origin = {0.0, 0.0, 0.0}
  amplitude = 1.0
  hwidth = 1.0
  return gauss_pulse(x, y, z, origin, amplitude, hwidth)
end

