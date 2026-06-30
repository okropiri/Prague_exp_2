#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MOUNT_POINT="${PADIWA_RAW_CAEN_MOUNT:-$SCRIPT_DIR/Raw_Caen_Data}"
REMOTE_USER="${PADIWA_USER:-padiwa}"
REMOTE_HOST="${PADIWA_HOST:-100.104.178.79}"
REMOTE_DIR="${PADIWA_REMOTE_DIR:-/data6/NCAL_Neutron_exp_REZ_cyclotron_v2}"
REMOTE_PATH="${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}"

usage() {
    cat <<EOF
Usage: $(basename "$0") [mount|unmount|remount|status]

Defaults:
  remote:      ${REMOTE_PATH}
  mount point: ${MOUNT_POINT}

Optional environment overrides:
  PADIWA_USER
  PADIWA_HOST
  PADIWA_REMOTE_DIR
  PADIWA_RAW_CAEN_MOUNT
EOF
}

require_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "Missing required command: $1" >&2
        exit 1
    fi
}

is_mounted() {
    mountpoint -q "$MOUNT_POINT"
}

ensure_mount_dir() {
    mkdir -p "$MOUNT_POINT"
}

ensure_empty_dir() {
    if [[ -n "$(find "$MOUNT_POINT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
        echo "Refusing to mount over a non-empty directory: $MOUNT_POINT" >&2
        exit 1
    fi
}

mount_remote() {
    require_cmd sshfs
    require_cmd mountpoint
    ensure_mount_dir

    if is_mounted; then
        echo "Already mounted: $MOUNT_POINT"
        return 0
    fi

    ensure_empty_dir
    sshfs "$REMOTE_PATH" "$MOUNT_POINT" -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3
    echo "Mounted $REMOTE_PATH at $MOUNT_POINT"
}

unmount_remote() {
    require_cmd mountpoint

    if ! is_mounted; then
        echo "Not mounted: $MOUNT_POINT"
        return 0
    fi

    if command -v fusermount >/dev/null 2>&1; then
        fusermount -u "$MOUNT_POINT"
    else
        umount "$MOUNT_POINT"
    fi

    echo "Unmounted $MOUNT_POINT"
}

show_status() {
    echo "Remote path: $REMOTE_PATH"
    echo "Mount point: $MOUNT_POINT"

    if is_mounted; then
        echo "Status: mounted"
        if command -v findmnt >/dev/null 2>&1; then
            findmnt --target "$MOUNT_POINT"
        fi
    else
        echo "Status: not mounted"
    fi
}

case "${1:-mount}" in
    mount)
        mount_remote
        ;;
    unmount)
        unmount_remote
        ;;
    remount)
        unmount_remote
        mount_remote
        ;;
    status)
        show_status
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac