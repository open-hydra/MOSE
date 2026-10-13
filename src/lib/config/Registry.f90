module MOSE_Input_Registry

  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  integer, parameter :: TYPE_INT=1, TYPE_REAL=2, TYPE_LOG=3, TYPE_STR=4

  public :: registry_t, Validate_Registry

  !--------------------------------------------------------
  ! Value container (typed pointers)
  !--------------------------------------------------------
  type :: param_value_t
    integer, pointer :: i => null()
    real(R8), pointer :: r => null()
    logical, pointer :: l => null()
    character(len=:), pointer :: s => null()

    integer, pointer :: iarr(:) => null()
    real(R8), pointer :: rarr(:) => null()
  end type


  !--------------------------------------------------------
  ! Parameter description
  !--------------------------------------------------------
  type :: param_t
    character(len=:), allocatable :: section
    character(len=:), allocatable :: name
    character(len=:), allocatable :: description
    character(len=:), allocatable :: default_str
    character(len=:), allocatable :: allowed

    logical :: required = .false.
    logical :: is_set   = .false.

    integer :: type_id = 0
    type(param_value_t) :: value
  end type


  !--------------------------------------------------------
  ! Registry container
  !--------------------------------------------------------
  type :: registry_t

    type(param_t), allocatable :: params(:)

    integer :: size = 0
    integer :: capacity = 0

  contains

    procedure :: reserve
    procedure :: add_int
    procedure :: add_real
    procedure :: add_logical
    procedure :: add_string

    procedure :: add_int_array
    procedure :: add_real_array

    procedure :: generate_markdown

    generic :: add => add_int, add_real, add_logical, add_string, &
        add_int_array, add_real_array

  end type


  type(registry_t), public :: reg


contains


  !========================================================
  ! Capacity management
  !========================================================

  subroutine reserve(this,newcap)

    implicit none

    class(registry_t), intent(inout) :: this
    integer, intent(in) :: newcap

    type(param_t), allocatable :: tmp(:)

    if (newcap <= this%capacity) return

    allocate(tmp(newcap))

    if (this%size > 0) then
        tmp(1:this%size) = this%params(1:this%size)
    end if

    call move_alloc(tmp,this%params)

    this%capacity = newcap

  end subroutine reserve


  !========================================================
  ! Internal: ensure space
  !========================================================

  subroutine ensure_space(this)

    implicit none

    class(registry_t), intent(inout) :: this

    if (this%size == this%capacity) then

      if (this%capacity == 0) then
        call this%reserve(8)
      else
        call this%reserve(2*this%capacity)
      end if

    end if

  end subroutine ensure_space


  !========================================================
  ! Add integer parameter
  !========================================================

  subroutine add_int(this,section,name,var,default,desc,allowed,required)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, target, intent(inout) :: var

    integer :: n

    call ensure_space(this)

    this%size = this%size + 1
    n = this%size

    this%params(n)%section = section
    this%params(n)%name = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed = allowed
    this%params(n)%required = required

    this%params(n)%type_id = TYPE_INT
    this%params(n)%value%i => var

    read(default,*) var

  end subroutine


  !========================================================
  ! Add real parameter
  !========================================================

  subroutine add_real(this,section,name,var,default,desc,allowed,required)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    real(R8), target, intent(inout) :: var

    integer :: n

    call ensure_space(this)

    this%size = this%size + 1
    n = this%size

    this%params(n)%section = section
    this%params(n)%name = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed = allowed
    this%params(n)%required = required

    this%params(n)%type_id = TYPE_REAL
    this%params(n)%value%r => var

    read(default,*) var

  end subroutine add_real


  !========================================================
  ! Add logical parameter
  !========================================================

  subroutine add_logical(this,section,name,var,default,desc,allowed,required)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    logical, target, intent(inout) :: var

    integer :: n

    call ensure_space(this)

    this%size = this%size + 1
    n = this%size

    this%params(n)%section = section
    this%params(n)%name = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed = allowed
    this%params(n)%required = required

    this%params(n)%type_id = TYPE_LOG
    this%params(n)%value%l => var

    read(default,*) var

  end subroutine add_logical


  !========================================================
  ! Add string parameter
  !========================================================

  subroutine add_string(this,section,name,var,default,desc,allowed,required)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    character(len=*), target, intent(inout) :: var

    integer :: n

    call ensure_space(this)

    this%size = this%size + 1
    n = this%size

    this%params(n)%section = section
    this%params(n)%name = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed = allowed
    this%params(n)%required = required

    this%params(n)%type_id = TYPE_STR
    this%params(n)%value%s => var

    var = default

  end subroutine add_string


  !========================================================
  ! Integer array
  !========================================================

  subroutine add_int_array(this,section,name,var,default,desc,allowed,required)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, target, intent(inout) :: var(:)

    integer :: n, defval

    call ensure_space(this)

    this%size = this%size + 1
    n = this%size

    this%params(n)%section = section
    this%params(n)%name = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed = allowed
    this%params(n)%required = required

    this%params(n)%type_id = TYPE_INT
    this%params(n)%value%iarr => var

    read(default,*) defval
    var(:) = defval

  end subroutine


  !========================================================
  ! Real array
  !========================================================

  subroutine add_real_array(this,section,name,var,default,desc,allowed,required)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    real(R8), target, intent(inout) :: var(:)

    integer :: n
    real(R8) :: defval

    call ensure_space(this)

    this%size = this%size + 1
    n = this%size

    this%params(n)%section = section
    this%params(n)%name = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed = allowed
    this%params(n)%required = required

    this%params(n)%type_id = TYPE_REAL
    this%params(n)%value%rarr => var

    read(default,*) defval
    var(:) = defval

  end subroutine add_real_array


  !========================================================
  ! Validation
  !========================================================

  !> Every rule of the column Allowed is checked: a bound (> x, >= x, < x, <= x) or an interval (in (a, b), with [ or ]
  !> for an end that belongs to it) for a number, on each value of a list of numbers; a list of names for a string;
  !> logical (or the list true, false) for a logical, which the reader of input.ini takes only as a logical. A rule of
  !> another form stops the run as well, so that no rule is skipped in silence.
  function Validate_Registry() result(out)

    implicit none

    character(len=1024) :: out
    integer :: i, k
    real(R8) :: val
    character(len=16) :: nth

    out=""

    do i=1,reg%size

      if (reg%params(i)%required .and. .not. reg%params(i)%is_set) then
        out = "[ERROR] Required parameter not set: "//trim(reg%params(i)%name)
        return
      end if

      if (reg%params(i)%allowed == "") cycle

      ! A list of numbers: the rule holds for each of its values
      if (associated(reg%params(i)%value%iarr)) then
        do k = 1, size(reg%params(i)%value%iarr)
          val = real(reg%params(i)%value%iarr(k), R8)
          call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
          if (out /= "") exit
        end do
      else if (associated(reg%params(i)%value%rarr)) then
        do k = 1, size(reg%params(i)%value%rarr)
          val = reg%params(i)%value%rarr(k)
          call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
          if (out /= "") exit
        end do
      end if
      if (associated(reg%params(i)%value%iarr) .or. associated(reg%params(i)%value%rarr)) then
        if (out /= "") then
          write(nth,'(I0)') k
          out = trim(out)//" (value "//trim(nth)//" of the list)"
          return
        end if
        cycle
      end if

      select case(reg%params(i)%type_id)

      case(TYPE_INT)

        if (.not. associated(reg%params(i)%value%i)) cycle

        val = real(reg%params(i)%value%i, R8)
        call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
        if (out /= "") return

      case(TYPE_REAL)

        if (.not. associated(reg%params(i)%value%r)) cycle

        val = reg%params(i)%value%r
        call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
        if (out /= "") return

      case(TYPE_STR)

        if (.not. associated(reg%params(i)%value%s)) cycle
        if (trim(reg%params(i)%value%s) == "") cycle
        call validate_string(reg%params(i)%name, reg%params(i)%value%s, reg%params(i)%allowed, out)
        if (out /= "") return

      case(TYPE_LOG)

        ! The reader takes a logical value only (a list-directed read: T or F, after an optional period, starts the value,
        ! and any other value stops the run): nothing is left to check, if the rule says no more than that
        if (.not. logical_rule(reg%params(i)%allowed)) then
          out = "[ERROR] "//trim(reg%params(i)%name)//": the registry does not check its rule "//trim(reg%params(i)%allowed)
          return
        end if

      case default
        cycle

      end select

    end do

  end function Validate_Registry


  !> The rule of a logical: logical, or the list of its two values (true, false; with or without the periods)
  logical function logical_rule(rule)

  implicit none

  character(*), intent(in) :: rule
  character(len=len(rule)) :: low
  character(len=:), allocatable :: r
  integer :: k

  low = lowercase(rule)
  r = ""
  do k = 1, len(low)
    if (low(k:k) /= " " .and. low(k:k) /= ".") r = r//low(k:k)
  end do
  logical_rule = r == "logical" .or. r == "true,false" .or. r == "false,true"

  end function logical_rule


  !> The rule of a number: a bound (> x, >= x, < x, <= x) or an interval, in (a, b), where ( and ) leave the end out and
  !> [ and ] take it in. A rule of another form, or whose numbers do not read, is an error of the registry itself.
  subroutine validate_numeric(name,val,rule,out)

  implicit none

  character(*), intent(in) :: name
  real(R8), intent(in) :: val
  character(*), intent(in) :: rule
  character(len=*), intent(inout) :: out

  real(R8) :: limit, lo, hi
  character(len=:), allocatable :: r
  integer :: c, n, ios1, ios2
  logical :: ok, known

  r = trim(adjustl(rule))
  n = len(r)
  known = .true.
  ok = .true.
  ios1 = 0
  ios2 = 0

  if (index(r,"in") == 1) then

    ! in (a, b) and the like: the first and last characters after "in" open and close the interval
    r = trim(adjustl(r(3:)))
    n = len(r)
    c = index(r,",")
    known = n >= 5 .and. c > 2 .and. c < n - 1
    if (known) known = scan(r(1:1),"([") == 1 .and. scan(r(n:n),")]") == 1
    if (known) then
      read(r(2:c-1),*,iostat=ios1) lo
      read(r(c+1:n-1),*,iostat=ios2) hi
      known = ios1 == 0 .and. ios2 == 0
    end if
    if (known) ok = merge(val >= lo, val > lo, r(1:1) == "[") .and. merge(val <= hi, val < hi, r(n:n) == "]")
    r = "in "//r

  else if (index(r,">=") == 1) then

    read(r(3:),*,iostat=ios1) limit
    ok = val >= limit

  else if (index(r,"<=") == 1) then

    read(r(3:),*,iostat=ios1) limit
    ok = val <= limit

  else if (index(r,">") == 1) then

    read(r(2:),*,iostat=ios1) limit
    ok = val > limit

  else if (index(r,"<") == 1) then

    read(r(2:),*,iostat=ios1) limit
    ok = val < limit

  else

    known = .false.

  end if

  if (.not. known .or. ios1 /= 0) then
    out = "[ERROR] "//trim(name)//": the registry does not check its rule "//trim(rule)
  else if (.not. ok) then
    out = "[ERROR] "//trim(name)//" must be "//r
  end if

  end subroutine validate_numeric


  subroutine validate_string(name,value,allowed,out)

  implicit none

  character(*), intent(in) :: name
  character(*), intent(in) :: value
  character(*), intent(in) :: allowed
  character(len=*), intent(inout) :: out
  character(len=:), allocatable :: token
  character(len=:), allocatable :: lower_value
  integer :: i, start, end
  logical :: is_sep

  ! Split allowed by comma and check each token
  start = 1
  do i = 1, len(allowed)+1
    is_sep = i > len(allowed)
    if (.not. is_sep) is_sep = (allowed(i:i) == ",")
    if (is_sep) then
        end = i - 1
        token = adjustl(allowed(start:end))
        ! [NOT ACTIIVE] Allow capital and lowercase values to compare them in a case-insensitive way
        ! token = lowercase(token)
        ! lower_value = lowercase(value)
        lower_value = value
        if (trim(token) == trim(lower_value)) then
          return
        end if
        start = i + 1
    end if
  end do

  ! If we get here, no match was found in the allowed list
  out = "[ERROR] "//trim(name)//" must be one of: "//trim(allowed)

  end subroutine validate_string


  !========================================================
  ! Markdown generator
  !========================================================

  subroutine generate_markdown(this,filename)

    implicit none

    class(registry_t), intent(in) :: this
    character(*), intent(in), optional :: filename

    integer :: i,unit
    character(len=:), allocatable :: fileout
    character(len=:), allocatable :: current_section
    logical :: new_section

    if (present(filename)) then
      fileout = filename
    else
      fileout = "registry.md"
    end if

    open(newunit=unit,file=fileout,status="replace")

    write(unit,'(A)') "# Input Parameters"
    write(unit,'(A)') ""

    do i=1,this%size

      if (.not. allocated(current_section)) then
        new_section = .true.
      else
        new_section = trim(this%params(i)%section) /= trim(current_section)
      end if

      if (new_section) then
        current_section = this%params(i)%section
        write(unit,'(A)') ""
        write(unit,'(A)') "## "//trim(current_section)
        write(unit,'(A)') ""
        write(unit,'(A)') "| Parameter | Default | Allowed | Required | Description |"
        write(unit,'(A)') "|-----------|---------|---------|----------|-------------|"
      end if

      write(unit,'(A)') "| "//trim(this%params(i)%name)// &
          " | "//trim(this%params(i)%default_str)// &
          " | "//trim(this%params(i)%allowed)// &
          " | "//merge("yes"," no",this%params(i)%required)// &
          " | "//trim(this%params(i)%description)//" |"

    end do

    close(unit)

  end subroutine generate_markdown


  function lowercase(str) result(out)

    implicit none

    character(*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i

    do i = 1, len(str)
        if (str(i:i) >= "A" .and. str(i:i) <= "Z") then
            out(i:i) = char(iachar(str(i:i)) + 32)
        else
            out(i:i) = str(i:i)
        end if
    end do

  end function lowercase


end module MOSE_Input_Registry