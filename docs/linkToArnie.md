To link a new package to Arnie, PR the following modifications to `bpps.py`:

Add to the `bpps` function:

```

def bpps(sequence, ...., myNewArgument=myDefaultValue)

   # New package wrapping
   if package=='myGreatPackage':
      return bpps_myGreatPackage_(sequence, myNewArgument=myNewArgument)
      
   # End new package wrapping


```

And define a new function in `bpps.py`:

```
      
def bpps_myGreatPackage_(sequence, myNewArgument=True):

   # wrap algorithm, return calculations
   
  return bp_array #bp_array is a symmetric matrix of probabilities p(i:j) as a numpy array
  
 ```
  

In `eternabench/package_options.py`: add a keyword to the dictionary `package_options` is fed to arnie to control the keywords used for each package call.

```

package_options = {
'vienna_2': {'package': 'vienna_2'},
'vienna_2_nodangles': {'package':'vienna_2', 'dangles': False},
'vienna_2_60C': {'package': 'vienna_2', 'T': 60},
'vienna_1': {'package': 'vienna_1'},

'myGreatPackage': {'package': 'myGreatPackage', 'myNewArgument':myArgumentValue}
```
