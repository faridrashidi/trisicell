Mutation Calling
----------------

The overview of the mutation calling pipeline that was used in Trisicell is provided below.

.. raw:: html
    
    <br>
    <img src="_static/images/gatk.png" align="center">
    <br>
    <br>

We provide a more detailed description of each of these steps.

Step 1:
~~~~~~~

.. code-block:: bash
    
    STAR \
        --runMode genomeGenerate \
        --genomeDir {outdir_indexing_1} \
        --genomeFastaFiles {ref} \
        --sjdbGTFfile {annot} \
        --sjdbOverhang {readlength} \
        --runThreadN {number_of_threads}
