class InsertionException(Exception):
    """
    Detected a site which contains an insertion
    """
    pass

class DeletionException(Exception):
    """
    Detected a site which contains a deletion
    """
    pass

class PolyallelicException(Exception):
    """
    Detected a site which contains more than two alleles
    """
    pass

class HeterozygousException(Exception):
    """
    Detected a site which is heterozygous
    """
    pass





class CoverageError(Exception):
    """
    Stated depth of coverage does not equal the length of the bases string
    """
    pass

class OutOfSyncPileupError(Exception):
    """
    The position values do not match across all pileups
    """
    pass


