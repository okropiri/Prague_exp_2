#define _GNU_SOURCE

#include <dlfcn.h>
#include <dirent.h>
#include <fcntl.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>

static const char *legacy_path = "/dev/usb/v1718_0";
static const char *vendor_id = "21e1";
static const char *product_id = "0000";

static int (*real_open_fn)(const char *pathname, int flags, ...);
static int (*real_open64_fn)(const char *pathname, int flags, ...);
static int (*real_openat_fn)(int dirfd, const char *pathname, int flags, ...);
static int (*real_openat64_fn)(int dirfd, const char *pathname, int flags, ...);

__attribute__((constructor))
static void init_syms(void) {
    if (!real_open_fn) {
        real_open_fn = dlsym(RTLD_NEXT, "open");
    }
    if (!real_open64_fn) {
        real_open64_fn = dlsym(RTLD_NEXT, "open64");
    }
    if (!real_openat_fn) {
        real_openat_fn = dlsym(RTLD_NEXT, "openat");
    }
    if (!real_openat64_fn) {
        real_openat64_fn = dlsym(RTLD_NEXT, "openat64");
    }
}

static int read_trimmed(const char *path, char *buffer, size_t size) {
    int fd = syscall(SYS_openat, AT_FDCWD, path, O_RDONLY);
    ssize_t count;

    if (fd < 0) {
        return -1;
    }

    count = read(fd, buffer, size - 1);
    close(fd);
    if (count <= 0) {
        return -1;
    }

    buffer[count] = '\0';
    buffer[strcspn(buffer, "\r\n")] = '\0';
    return 0;
}

static int find_caen_usb(char *out, size_t out_size) {
    DIR *dir = opendir("/sys/bus/usb/devices");
    if (!dir) {
        return -1;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        char path[PATH_MAX];
        char idv[32];
        char idp[32];
        char bus[32];
        char dev[32];

        if (entry->d_name[0] == '.') {
            continue;
        }

        snprintf(path, sizeof(path), "/sys/bus/usb/devices/%s/idVendor", entry->d_name);
        if (read_trimmed(path, idv, sizeof(idv)) != 0) {
            continue;
        }

        snprintf(path, sizeof(path), "/sys/bus/usb/devices/%s/idProduct", entry->d_name);
        if (read_trimmed(path, idp, sizeof(idp)) != 0) {
            continue;
        }

        if (strcmp(idv, vendor_id) != 0 || strcmp(idp, product_id) != 0) {
            continue;
        }

        snprintf(path, sizeof(path), "/sys/bus/usb/devices/%s/busnum", entry->d_name);
        if (read_trimmed(path, bus, sizeof(bus)) != 0) {
            continue;
        }

        snprintf(path, sizeof(path), "/sys/bus/usb/devices/%s/devnum", entry->d_name);
        if (read_trimmed(path, dev, sizeof(dev)) != 0) {
            continue;
        }

        snprintf(out, out_size, "/dev/bus/usb/%03d/%03d", atoi(bus), atoi(dev));
        closedir(dir);
        return 0;
    }

    closedir(dir);
    return -1;
}

static const char *redirect_path(const char *pathname, char *buffer, size_t size) {
    if (!pathname) {
        return pathname;
    }
    if (strcmp(pathname, legacy_path) != 0) {
        return pathname;
    }
    if (find_caen_usb(buffer, size) == 0) {
        return buffer;
    }
    return pathname;
}

static int do_open(int (*fn)(const char *, int, ...), const char *pathname, int flags, va_list args) {
    char replacement[PATH_MAX];
    const char *target;

    init_syms();
    target = redirect_path(pathname, replacement, sizeof(replacement));

    if (flags & O_CREAT) {
        mode_t mode = va_arg(args, mode_t);
        return fn(target, flags, mode);
    }
    return fn(target, flags);
}

static int do_openat(int (*fn)(int, const char *, int, ...), int dirfd, const char *pathname, int flags, va_list args) {
    char replacement[PATH_MAX];
    const char *target;

    init_syms();
    target = redirect_path(pathname, replacement, sizeof(replacement));

    if (flags & O_CREAT) {
        mode_t mode = va_arg(args, mode_t);
        return fn(dirfd, target, flags, mode);
    }
    return fn(dirfd, target, flags);
}

int open(const char *pathname, int flags, ...) {
    int rc;
    va_list args;

    va_start(args, flags);
    rc = do_open(real_open_fn, pathname, flags, args);
    va_end(args);
    return rc;
}

int open64(const char *pathname, int flags, ...) {
    int rc;
    va_list args;

    va_start(args, flags);
    rc = do_open(real_open64_fn, pathname, flags, args);
    va_end(args);
    return rc;
}

int openat(int dirfd, const char *pathname, int flags, ...) {
    int rc;
    va_list args;

    va_start(args, flags);
    rc = do_openat(real_openat_fn, dirfd, pathname, flags, args);
    va_end(args);
    return rc;
}

int openat64(int dirfd, const char *pathname, int flags, ...) {
    int rc;
    va_list args;

    va_start(args, flags);
    rc = do_openat(real_openat64_fn, dirfd, pathname, flags, args);
    va_end(args);
    return rc;
}