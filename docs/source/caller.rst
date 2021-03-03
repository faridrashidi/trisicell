Mutation Calling
----------------

Our mutation calling procedures has the following steps:

.. raw:: html
    
    <br>
    <img src="_static/images/gatk.png" align="center">
    <br>
    <br>

In the first step:

.. code-block:: bash
    
    STAR \
        --runMode genomeGenerate \
        --genomeDir {outdir_indexing_1} \
        --genomeFastaFiles {ref} \
        --sjdbGTFfile {annot} \
        --sjdbOverhang {readlength} \
        --runThreadN {number_of_threads}
