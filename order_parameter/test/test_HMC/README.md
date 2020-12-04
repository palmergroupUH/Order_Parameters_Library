## Validating the HMC order parameter 

### Test trajectory description:

* 512 mW coarse-grained model of water

* T = 200 K and P = 1 bar

* 50 ns production run

* Total of 50 configuratinos sampled with 1 ns interval 

* Spontaneous crystallization 

### Instructions:

* If the pytest is installed

```
pytest 
``` 

* otherwise, 

```
python test_order_parameter.py  
```
