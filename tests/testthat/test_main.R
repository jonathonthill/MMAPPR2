context("Main and helper functions")

vepFlags <- readRDS('test_data/objects/vep_flags.RDS')
param <- MmapprParam(refGenome=new("GmapGenome"),
                     wtFiles='test_data/bam_files/zy14_dummy.bam',
                     mutFiles='test_data/bam_files/zy14_dummy.bam',
                     species='danio_rerio',
                     vepFlags=vepFlags,
                     refFasta='test_data/dummy.fasta')


test_that("MmapprParam has default VEPFlags when one isn't provided", {
    # put fake vep in path since VEPFlags requires it
    dir.create('/tmp/ensembl-vep')
    system('touch /tmp/ensembl-vep/vep; chmod 777 /tmp/ensembl-vep/vep')
    originalPath <- Sys.getenv('PATH')
    Sys.setenv('PATH'=paste0(originalPath, ':', '/tmp/ensembl-vep'))

    mp <- MmapprParam(refGenome=new("GmapGenome"),
                      wtFiles='test_data/bam_files/zy14_dummy.bam',
                      mutFiles='test_data/bam_files/zy14_dummy.bam',
                      species='danio_rerio',
                      refFasta='test_data/dummy.fasta')
    expect_false(is.null(mp@vepFlags))
    expect_true(ensemblVEP::flags(vepFlags(mp))$database == FALSE)
    expect_true(ensemblVEP::flags(vepFlags(mp))$species == 'danio_rerio')

    # cleanup
    Sys.setenv('PATH'=originalPath)
    unlink('/tmp/ensembl-vep', recursive=TRUE)
})


test_that("MmapprParam takes only real files; handles character or objects", {
    fn_wt <- 'test_data/bam_files/zy14_dummy.bam'
    fn_mut <- 'test_data/bam_files/zy14_dummy.bam'

    # non-existing filename shouldn't work
    expect_error(MmapprParam(refGenome=new("GmapGenome"), wtFiles="wt", mutFiles='mut',
                             species='danio_rerio', vepFlags=vepFlags,
                             refFasta='test_data/dummy.fasta'))

    # existing filename, just character, should work
    param <- MmapprParam(refGenome=new("GmapGenome"), wtFiles=fn_wt, mutFiles=fn_mut,
                         species='danio_rerio', vepFlags=vepFlags,
                         refFasta='test_data/dummy.fasta')
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    # character with length>1
    param <- MmapprParam(refGenome=new("GmapGenome"), wtFiles=c(fn_wt, fn_wt), mutFiles=fn_mut,
                         species='danio_rerio', vepFlags=vepFlags,
                         refFasta='test_data/dummy.fasta')
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    bf <- Rsamtools::BamFile(fn_mut)
    param <- MmapprParam(refGenome=new("GmapGenome"), wtFiles=bf, mutFiles=bf,
                         species='danio_rerio', vepFlags=vepFlags,
                         refFasta='test_data/dummy.fasta')
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    bfl <- Rsamtools::BamFileList(fn_mut)
    param <- MmapprParam(refGenome=new("GmapGenome"), wtFiles=bfl, mutFiles=bfl,
                         species='danio_rerio', vepFlags=vepFlags,
                         refFasta='test_data/dummy.fasta')
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    expect_error(MmapprParam(refGenome=new("GmapGenome"), wtFiles=bfl, mutFiles=bfl,
                             species='danio_rerio', vepFlags=vepFlags,
                             refFasta='fake.fasta'))

})

test_that('file setters for param can change format automatically', {
    fn <- 'test_data/bam_files/zy14_dummy.bam'
    param <- MmapprParam(refGenome=new("GmapGenome"), wtFiles=fn, mutFiles=fn,
                         species='danio_rerio', vepFlags=vepFlags,
                         refFasta='test_data/dummy.fasta')

    wtFiles(param) <- c(fn, fn)
    mutFiles(param) <- c(fn, fn)
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    bf <- Rsamtools::BamFile(fn)
    wtFiles(param) <- c(bf, bf)
    mutFiles(param) <- c(bf, bf)
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")
})

# test_that('log writes to default folder', {
#     skip_if_not_installed('mockery')
#     newParam <- param
#     unlink('/tmp/m2', recursive=TRUE)
#     dir.create('/tmp/m2')
#     newParam@outputFolder <- '/tmp/m2'
#     mockery::stub(mmappr,
#                   '.prepareOutputFolder',
#                   new('MmapprData', param = newParam))
#     mockery::stub(mmappr, 'tryCatch', function(...) stop('expected failure'))
#     expect_error(mmappr(param)) # , 'expected failure'
#     expect_true(file.exists('/tmp/m2/mmappr2.log'))
#     unlink('/tmp/m2', recursive=TRUE)
# })

test_that('checkDep finds program iff in path', {
    tryCatch({
        tempProgramLocation <- tempdir()
        originalPath = Sys.getenv('PATH')
        Sys.setenv(PATH=tempProgramLocation)
        file.create(file.path(tempProgramLocation, 'program'))
        Sys.chmod(file.path(tempProgramLocation, 'program'))
        expect_true(.checkDep('program'))
        expect_error(.checkDep('rubbish'))
        Sys.setenv(PATH = 'nonexistentFolder')
        expect_error(.checkDep('program'))
    }, finally = {
        Sys.setenv(PATH=originalPath)
    })

})
