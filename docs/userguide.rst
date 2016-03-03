User manual
==================================

   
The mrchem input file
---------------------

The input file is organized in sections and keywords that can be of different
type 

.. code-block:: c
    
     Section {
       keyword_1 = 1
       keyword_2 = 3.14
       keyword_3 = [1, 2, 3]
       keyword_4 = "foo"
       keyword_5 = True
    }

.. highlight:: c
    
     Section {
       keyword_1 = 1
       keyword_2 = 3.14
       keyword_3 = [1, 2, 3]
       keyword_4 = "foo"
       keyword_5 = True
    }


The main input section contain two important keywords that specify the
polynomial order of the multiwavelet basis set, and the relative precision that
will be guaranteed in the calculation. The main input section is not specified
by name, just write the keywords directly, e.g


