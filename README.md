# PeakPickingTools

PeakPickingTools is a MATLAB toolbox to exctract modal parameters from frequency response functions, using a robust peak-picking identification technique.

It is an implementation of the method described in [[1]](#references).

<img src="logo.png" alt="PeakPickingTools logo" width="300"/>


## Installation and usage

Download or clone the repository. Folder containing the toolbox has to be set as current folder in MATLAB, or added to path.

Run .m file.

```
>> PeakPickingTools
```

Input data (input impedance or any FRF) must be provided as a text file containing three columns : frequency, real part, imaginary part. The text file must have .dat or .txt extension.


## Miscellaneous

* PeakPickingTools is not yet documented.
* Bug reports and general feedback are welcome ([frederic.ablitzer@univ-lemans.fr](mailto:frederic.ablitzer@univ-lemans.fr?subject=PeakPickingTools)). You may consider sending data as attached file.
* If used in research, please cite the original research paper for which the code was developed for [[1]](#references).
* PeakPickingTools has been developped and tested in MATLAB version R2020a. Compatibility with other versions is not guaranteed.

## Version history

* 1.0
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.



## References

[1] Ablitzer F. 2021. Peak-picking identification technique for modal expansion of input impedance of brass instruments. Acta Acustica, 5, 53. [https://doi.org/10.1051/aacus/2021046](https://doi.org/10.1051/aacus/2021046)
