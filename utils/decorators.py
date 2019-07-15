from functools import wraps


def exception_catcher(func):
    """
    basically just used to catch and ignore exceptions when running in bulk.
    This prevents the whole program stopping when just one molecule is 'broken'.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):

        try:
            return func(*args, **kwargs)

        except Exception as exc:
            if len(args) >= 1 and hasattr(args[0], 'molecule'):
                if hasattr(args[0].molecule, 'bulk_run'):
                    if not args[0].molecule.bulk_run:
                        raise
                else:
                    raise

            print(exc)

    return wrapper
