# class-matrix_class-sole-summer-practice-project

![Matrix](/pics/Matrix.gif)

![Static Badge](https://img.shields.io/badge/C%2B%2B-purple?style=flat&logo=C%2B%2B&label=Language)
![Static Badge](https://img.shields.io/badge/VS_Code-purple?style=flat&logo=VS%20Code&label=IDE&labelColor=darkblue&color=darkgreen)
![Static Badge](https://img.shields.io/badge/MinGW-purple?style=flat&logo=MinGW&label=env&labelColor=darkred&color=teal)


## Description
This project is an implementation of modules for working with matrices, as well as solutions to <u>____systems of linear equations____</u> that are implemented through <u>__matrices__</u> (one can say, the implementation of an online calculator for solving <u>***SOLE***</u> and operations on matrices). The implementation is presented in ***C++***.

## Features
- <u>__Execution of C++__</u> code in a convenient way
- **Integrating** the code into the application
- **Implementation of matrices** of complex numbers and <u>__a large number of methods__</u> for working with matrices and finding solutions to various __SOLE__


## Documentation
The implementation is presented in the directory [Scripts](https://github.com/BorDch/class-matrix_class-sole-summer-practice-project/tree/main/Scripts) in the __modules and source files__.

* __Modules__:
    * ___matrix.hpp___: header which includes ***Matrix class declaration*** (constructors, methods);

    * ___sole.hpp___: header which includes declaration ***of methods for solving systems***. 
* __Source__:
    * ___matrix.cpp___: definition and realization of __Matrix__ class from `matrix.hpp` 

    * ___matrix_main_cpp___:
    launching programm with linking from `matrix.hpp` and `matrix.cpp` with examples of launching.
    
    * ___sole.cpp___:
    definition and realization of __SOLE__ class from `sole.hpp`
    
    * ___sole_main.cpp___:
    launching programm with linking from `sole.hpp` and `sole.cpp` with examples of launching.

## Installation
1. __Clone repository__:
    ```sh
    git clone https://github.com/BorDch/class-matrix_class-sole-summer-practice-project.git
    ```
2. __Go to the project directory__:
    ```sh
    cd class-matrix_class-sole-summer-practice-project
    ```
3. __Project Building__:
    ```sh
    mkdir build && cd build
    cmake ...
    make
    ```
4. __Execute__:
    ```sh
    ./exe_file_name
    ```

## Developers
- [Boris Cherkasov](https://github.com/BorDch)

## Useful Links
- [Matrices](https://en.wikipedia.org/wiki/Matrix_(mathematics))
- [SOLE (System of linear equations)](https://en.wikipedia.org/wiki/System_of_linear_equations)
- [Numerical methods for solving SOLE](https://docs.yandex.by/docs/view?tm=1738581625&tld=by&lang=en&name=nm-ang.pdf&text=All%20Methods%20for%20solving%20system%20of%20linear%20equations&url=https%3A%2F%2Fcw.fel.cvut.cz%2Fold%2F_media%2Fcourses%2Fa4b01num%2Fen%2Flectures%2Fnm-ang.pdf&lr=157&mime=pdf&l10n=ru&sign=9172ba8c3d8adffeac5894020512a87e&keyno=0&nosw=1&serpParams=tm%3D1738581625%26tld%3Dby%26lang%3Den%26name%3Dnm-ang.pdf%26text%3DAll%2BMethods%2Bfor%2Bsolving%2Bsystem%2Bof%2Blinear%2Bequations%26url%3Dhttps%253A%2F%2Fcw.fel.cvut.cz%2Fold%2F_media%2Fcourses%2Fa4b01num%2Fen%2Flectures%2Fnm-ang.pdf%26lr%3D157%26mime%3Dpdf%26l10n%3Dru%26sign%3D9172ba8c3d8adffeac5894020512a87e%26keyno%3D0%26nosw%3D1)

## Contributing
We will be glad to receive your <u>__PR (pull requests)</u>__! :) 

To make changes: 

1. __Fork repository__
2. __Create a branch__ (`git checkout -b  feature-branch`)
3. __Make changes and save them__ (`git commit -m "Add new features"`)
4. __Push them__ (`git push origin feature-branch`)
5. __Open PR__.

## License
The project is distributed under the `MIT` license.
