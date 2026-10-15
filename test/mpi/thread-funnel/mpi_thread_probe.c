/*
 * MPI thread-level probe for the MOSE MPI tests (test/mpi/thread-funnel).
 *
 * Preloaded into the solver (LD_PRELOAD), it wraps the Fortran MPI entry points
 * MOSE calls, remembers the thread that called MPI_INIT_THREAD (the "main
 * thread" of the MPI standard) and counts, per routine, how many calls come
 * from any other thread. Under MPI_THREAD_FUNNELED that count must be zero
 * (MPI 4.1, Sec. 11.2.1 and 11.6.2). It also keeps the largest tag of the
 * point-to-point calls (-1 if none): valid tags run from 0 to MPI_TAG_UB, and
 * the standard guarantees only MPI_TAG_UB >= 32767 (Sec. 3.2.3). The probe
 * makes no MPI call of its own; it identifies threads with pthread_self() and
 * forwards every call to the real routine found with dlsym(RTLD_NEXT).
 *
 * Just before the real MPI_FINALIZE each process prints one line to stderr:
 *   MPI-THREAD-PROBE rank=<r> required=<q> provided=<p> calls=<n> other=<m> maxtag=<t> [name:calls/other ...]
 * <r> comes from the launcher's environment (PMI_RANK, OMPI_COMM_WORLD_RANK,
 * PMIX_RANK), "?" if none is set.
 *
 * The wrapped names are the Fortran bindings with the usual lower-case +
 * trailing-underscore mangling of gfortran and ifx, and the same names without
 * the underscore, which a solver compiled with -fno-underscoring (the GNU
 * flags of MOSE) calls; every argument of these bindings is passed by
 * reference, so each wrapper forwards exactly the pointers of the binding it
 * replaces.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef void *P;

enum {
  I_ABORT, I_ALLREDUCE, I_BARRIER, I_BCAST, I_COMM_RANK, I_COMM_SIZE, I_FINALIZE,
  I_GATHER, I_INITIALIZED, I_INIT_THREAD, I_IRECV, I_ISEND, I_IS_THREAD_MAIN,
  I_QUERY_THREAD, I_RECV_INIT, I_REDUCE, I_REQUEST_FREE, I_SEND_INIT, I_STARTALL,
  I_WAITALL, I_WTIME, N_SYM
};

static const char *names[N_SYM] = {
  "mpi_abort_", "mpi_allreduce_", "mpi_barrier_", "mpi_bcast_", "mpi_comm_rank_",
  "mpi_comm_size_", "mpi_finalize_", "mpi_gather_", "mpi_initialized_",
  "mpi_init_thread_", "mpi_irecv_", "mpi_isend_", "mpi_is_thread_main_",
  "mpi_query_thread_", "mpi_recv_init_", "mpi_reduce_", "mpi_request_free_",
  "mpi_send_init_", "mpi_startall_", "mpi_waitall_", "mpi_wtime_"
};

static _Atomic(void *) reals[N_SYM];
static atomic_long calls[N_SYM], other[N_SYM];
static pthread_t main_thread;
static atomic_int main_known = 0;
static int required = -1, provided = -1;
static atomic_int max_tag = -1;

/* The real routine, looked up on first use; atomic because several threads
 * may make their first call at the same time (they store the same pointer). */
static void *real(int i)
{
  void *p = atomic_load(&reals[i]);
  if (!p) {
    p = dlsym(RTLD_NEXT, names[i]);
    if (!p) {
      fprintf(stderr, "MPI-THREAD-PROBE: symbol %s not found\n", names[i]);
      abort();
    }
    atomic_store(&reals[i], p);
  }
  return p;
}

static void count(int i)
{
  atomic_fetch_add(&calls[i], 1);
  if (atomic_load(&main_known) && !pthread_equal(pthread_self(), main_thread))
    atomic_fetch_add(&other[i], 1);
}

/* tag: the fifth argument of the point-to-point routines wrapped below */
static void note_tag(P tag)
{
  int t = *(int *)tag, m = atomic_load(&max_tag);
  while (t > m && !atomic_compare_exchange_weak(&max_tag, &m, t))
    ;
}

/* The line is built whole and written with a single write(): the launcher
 * merges the stderr of all ranks, and the pieces of separate writes of two
 * ranks can interleave into one line. A write of fewer than PIPE_BUF bytes
 * to a pipe is atomic (POSIX); the line is shorter than 1 kB. */
static void report(void)
{
  const char *r = getenv("PMI_RANK");
  if (!r) r = getenv("OMPI_COMM_WORLD_RANK");
  if (!r) r = getenv("PMIX_RANK");
  if (!r) r = "?";
  long n = 0, m = 0;
  for (int i = 0; i < N_SYM; i++) { n += atomic_load(&calls[i]); m += atomic_load(&other[i]); }
  char line[2048];
  size_t len = 0;
  len += snprintf(line + len, sizeof line - len,
                  "MPI-THREAD-PROBE rank=%s required=%d provided=%d calls=%ld other=%ld maxtag=%d",
                  r, required, provided, n, m, atomic_load(&max_tag));
  for (int i = 0; i < N_SYM && len < sizeof line; i++) {
    long c = atomic_load(&calls[i]);
    if (c) len += snprintf(line + len, sizeof line - len, " %s:%ld/%ld", names[i], c, atomic_load(&other[i]));
  }
  if (len > sizeof line - 1) len = sizeof line - 1;
  line[len++] = '\n';
  ssize_t w = write(STDERR_FILENO, line, len);
  (void)w;
}

#define FWD2(f, i) void f(P a, P b) { count(i); ((void (*)(P, P))real(i))(a, b); }
#define FWD3(f, i) void f(P a, P b, P c) { count(i); ((void (*)(P, P, P))real(i))(a, b, c); }
#define FWD4(f, i) void f(P a, P b, P c, P d) { count(i); ((void (*)(P, P, P, P))real(i))(a, b, c, d); }
#define FWD6(f, i) void f(P a, P b, P c, P d, P e, P g) { \
    count(i); ((void (*)(P, P, P, P, P, P))real(i))(a, b, c, d, e, g); }
#define FWD7(f, i) void f(P a, P b, P c, P d, P e, P g, P h) { \
    count(i); ((void (*)(P, P, P, P, P, P, P))real(i))(a, b, c, d, e, g, h); }
#define FWD8(f, i) void f(P a, P b, P c, P d, P e, P g, P h, P k) { \
    count(i); ((void (*)(P, P, P, P, P, P, P, P))real(i))(a, b, c, d, e, g, h, k); }
#define FWD8T(f, i) void f(P a, P b, P c, P d, P e, P g, P h, P k) { \
    count(i); note_tag(e); ((void (*)(P, P, P, P, P, P, P, P))real(i))(a, b, c, d, e, g, h, k); }
#define FWD9(f, i) void f(P a, P b, P c, P d, P e, P g, P h, P k, P l) { \
    count(i); ((void (*)(P, P, P, P, P, P, P, P, P))real(i))(a, b, c, d, e, g, h, k, l); }

/* Argument counts are those of the Fortran bindings, IERROR included. */
FWD3(mpi_abort_,          I_ABORT)          /* comm, errorcode, ierror                          */
FWD7(mpi_allreduce_,      I_ALLREDUCE)      /* sendbuf, recvbuf, count, type, op, comm, ierror  */
FWD2(mpi_barrier_,        I_BARRIER)        /* comm, ierror                                     */
FWD6(mpi_bcast_,          I_BCAST)          /* buf, count, type, root, comm, ierror             */
FWD3(mpi_comm_rank_,      I_COMM_RANK)      /* comm, rank, ierror                               */
FWD3(mpi_comm_size_,      I_COMM_SIZE)      /* comm, size, ierror                               */
FWD9(mpi_gather_,         I_GATHER)         /* sbuf, scount, stype, rbuf, rcount, rtype, root, comm, ierror */
FWD2(mpi_initialized_,    I_INITIALIZED)    /* flag, ierror                                     */
FWD8T(mpi_irecv_,         I_IRECV)          /* buf, count, type, source, tag, comm, req, ierror */
FWD8T(mpi_isend_,         I_ISEND)          /* buf, count, type, dest, tag, comm, req, ierror   */
FWD2(mpi_is_thread_main_, I_IS_THREAD_MAIN) /* flag, ierror                                     */
FWD2(mpi_query_thread_,   I_QUERY_THREAD)   /* provided, ierror                                 */
FWD8T(mpi_recv_init_,     I_RECV_INIT)      /* buf, count, type, source, tag, comm, req, ierror */
FWD8(mpi_reduce_,         I_REDUCE)         /* sendbuf, recvbuf, count, type, op, root, comm, ierror */
FWD2(mpi_request_free_,   I_REQUEST_FREE)   /* request, ierror                                  */
FWD8T(mpi_send_init_,     I_SEND_INIT)      /* buf, count, type, dest, tag, comm, req, ierror   */
FWD3(mpi_startall_,       I_STARTALL)       /* count, requests, ierror                          */
FWD4(mpi_waitall_,        I_WAITALL)        /* count, requests, statuses, ierror                */

double mpi_wtime_(void)
{
  count(I_WTIME);
  return ((double (*)(void))real(I_WTIME))();
}

void mpi_init_thread_(P req, P prov, P ierr)
{
  count(I_INIT_THREAD);
  ((void (*)(P, P, P))real(I_INIT_THREAD))(req, prov, ierr);
  required = *(int *)req;
  provided = *(int *)prov;
  main_thread = pthread_self();
  atomic_store(&main_known, 1);
}

void mpi_finalize_(P ierr)
{
  count(I_FINALIZE);
  report();
  ((void (*)(P))real(I_FINALIZE))(ierr);
}

/* The names without the trailing underscore (-fno-underscoring): the same wrappers, which forward to the routines of
 * the names with it (the MPI library exports both for the same routine). */
#define NOUNDER(f) extern __typeof__(f##_) f __attribute__((alias(#f "_")));
NOUNDER(mpi_abort)
NOUNDER(mpi_allreduce)
NOUNDER(mpi_barrier)
NOUNDER(mpi_bcast)
NOUNDER(mpi_comm_rank)
NOUNDER(mpi_comm_size)
NOUNDER(mpi_finalize)
NOUNDER(mpi_gather)
NOUNDER(mpi_initialized)
NOUNDER(mpi_init_thread)
NOUNDER(mpi_irecv)
NOUNDER(mpi_isend)
NOUNDER(mpi_is_thread_main)
NOUNDER(mpi_query_thread)
NOUNDER(mpi_recv_init)
NOUNDER(mpi_reduce)
NOUNDER(mpi_request_free)
NOUNDER(mpi_send_init)
NOUNDER(mpi_startall)
NOUNDER(mpi_waitall)
NOUNDER(mpi_wtime)
