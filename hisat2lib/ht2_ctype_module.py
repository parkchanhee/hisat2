#!/usr/bin/env python3

import ctypes


class HT2_OPTION(ctypes.Structure):
    """
    struct ht2_options {
        int offRate;

        int useMm;
        int useShmem;
        int mmSweep;
        int noRefNames;
        int noSplicedAlignment;
        int gVerbose;
        int startVerbose;
        int sanityCheck;

        int useHaplotype;
    };
    """
    _fields_ = [("offRate", ctypes.c_int),
                ("useMm", ctypes.c_int),
                ("useShmem", ctypes.c_int),
                ("mmSweep", ctypes.c_int),
                ("noRefNames", ctypes.c_int),
                ("noSplicedAlignment", ctypes.c_int),
                ("gVerbose", ctypes.c_int),
                ("startVerbose", ctypes.c_int),
                ("sanityCheck", ctypes.c_int),
                ("useHaplotype", ctypes.c_int)
                ]


class HT2_INDEX_GETREFNAMES_RESULT(ctypes.Structure):
    """
    struct ht2_index_getrefnames_result {
        int count;
        char* names[0];
    };
    """
    _fields_ = [("count", ctypes.c_int),
                ("names", ctypes.c_char_p*0)
                ]


ht2lib = ctypes.CDLL('../cmake-build-debug/libhisat2lib.so')
print(ht2lib)


libc = ctypes.CDLL('libc.so.6')

print(libc)
print(libc.time(None))
libc.printf(b'Hello World ctypes\n')
libc.printf(b'Hello World ctypes\n')

a = ctypes.c_bool(True)

print(a.value)

print(ctypes.c_int(-1).value)
print(ctypes.c_uint(-1).value)


s = "Hello World"

c_s = ctypes.c_wchar_p(s)
print(c_s.value)
c_s.value = "Other world"
print(c_s.value)
print(s)


printf = libc.printf
printf(b"Hello %f\n", ctypes.c_double(1.0))


printf.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_double]
ret=printf(b"Hello New World %s %d %f\n", b"Hello", -1, 10.0)
print('Return from printf:', ret)

print(ord("a"))


i = ctypes.c_int()

print(i)
i.value = 0xfffffffff
print(i.value)
f = ctypes.c_float()
s = ctypes.create_string_buffer(b'\000' * 32)
print(repr(s.value))

libc.sscanf(b"1 3.14 Hello", b"%d %f %s", ctypes.byref(i), ctypes.byref(f), s)

print(i.value, f.value, s.value)
print(repr(s.value))


# exit(0)
# return types

ht2option = HT2_OPTION()
print(ht2option.offRate)
print(ht2option.useMm)
ht2option.gVerbose = 1
ht2option.startVerbose = 1
print(ht2option.gVerbose)

ht2lib.ht2_init.restype = ctypes.c_void_p
handle = ht2lib.ht2_init(b'/home/parkch/work/00-data/indexes-grch38-snp144/HISAT2_22/22_rep', ctypes.pointer(ht2option))

print('{:x}'.format(handle))
print(type(handle))

ht2lib.ht2_index_getrefnamebyid.restype = ctypes.c_char_p
ht2lib.ht2_index_getrefnamebyid.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
result = ht2lib.ht2_index_getrefnamebyid(handle, ctypes.c_uint32(0))

print('Result from getrefnamebyid', result)

resultptr = ctypes.c_void_p()
print('Before calling', resultptr.value)
ht2lib.ht2_index_getrefnames.restype = ctypes.c_int
ht2lib.ht2_index_getrefnames.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
result = ht2lib.ht2_index_getrefnames(handle, ctypes.pointer(resultptr))
print('Result from getrefnames', result)
print('Value of resultptr {:x}'.format(resultptr.value))

# resultptr.value is the address of HT2_INDEX_REFNAMES_RESULT structure
namesresult = ctypes.POINTER(HT2_INDEX_GETREFNAMES_RESULT)
namesresult.value = ctypes.cast(resultptr, ctypes.POINTER(HT2_INDEX_GETREFNAMES_RESULT))
print(namesresult.contents.count)
for a in namesresult.contents.names:
    print(a)


ht2lib.ht2_close(ctypes.c_void_p(handle))
