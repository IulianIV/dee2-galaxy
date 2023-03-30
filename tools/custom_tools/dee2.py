import sys
from dee2_conn import DEE2


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    selector, data_set, species, out_file, load_funcs, load_func_zip = '', '', '', '', '', ''
    load_func = ''

    try:
        selector = sys.argv[1]
        data_set = sys.argv[2]
        species = sys.argv[3]
        out_file = sys.argv[4]
        load_funcs = sys.argv[5]
        load_func_zip = sys.argv[6]
        # use this variable to grab a getDEE2 function selected from galaxy
        # function = sys.argv[4]

    except Exception:
        stop_err("Usage: python dee2.py data_set_type species")

    dee2 = DEE2(supress_r_warnings=False)

    dee2.data_set = data_set
    dee2.species = species

    if out_file != 'false':
        dee2.outfile = out_file

    if load_funcs != 'false':
        try:
            load_func = getattr(dee2, load_funcs)(load_func_zip)
        except AttributeError:
            return print("Called load function does not exist.")

    get_dee2 = dee2.queryDEE2()
    print(f'final results:\n{get_dee2}')


if __name__ == '__main__':
    main()
