module MOSE_Input_Registry

  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  integer, parameter :: TYPE_INT=1, TYPE_REAL=2, TYPE_LOG=3, TYPE_STR=4

  ! What a reading of input.ini during the run does with a key: the class that every reg%add gives (its argument during_run,
  ! the column During the run of docs/user/registry.md)
  integer, parameter, public :: KEY_TAKEN = 1      ! Every reading takes the value written
  integer, parameter, public :: KEY_FIXED = 2      ! Fixed for the run: the set-up makes from it a choice that a reading cannot
                                                   ! make again; a value written other than the one the set-up read is an error
                                                   ! of the reading (Fixed_Key_Changed)
  integer, parameter, public :: KEY_FIXED_OWN = 3  ! Fixed for the run, and checked by its model, which compares the value read
                                                   ! again by its meaning and names the error (the rotating frame)
  integer, parameter, public :: KEY_START = 4      ! For the start of the run only: a reading leaves it as it was

  public :: registry_t, Validate_Registry, Save_Values, Restore_Values, Key_Is_Set
  public :: Save_Start_Values, Fixed_Key_Changed, Keep_Fixed_Values

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
    logical :: unread   = .false.               ! In input.ini with a value that does not read as its type (Load_Ini)
    character(len=:), allocatable :: text       ! The value in input.ini as written, when the key is there (Load_Ini)

    integer :: type_id = 0
    integer :: during_run = 0                   ! KEY_TAKEN, KEY_FIXED, KEY_FIXED_OWN or KEY_START
    type(param_value_t) :: value
  end type


  !--------------------------------------------------------
  ! Registry container
  !--------------------------------------------------------
  type :: registry_t

    type(param_t), allocatable :: params(:)

    integer :: size = 0
    integer :: capacity = 0

    ! An error of input.ini outside the keys of the registry, which Validate_Registry gives first (Scan_Ini: levels)
    character(len=:), allocatable :: load_error

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


  !--------------------------------------------------------
  ! Copy of the values of the registry (Save_Values, Restore_Values; Save_Start_Values)
  !--------------------------------------------------------
  type :: saved_value_t
    integer  :: i = 0
    real(R8) :: r = 0.0_R8
    logical  :: l = .false.
    character(len=:), allocatable :: s
    integer,  allocatable :: iarr(:)
    real(R8), allocatable :: rarr(:)
    logical  :: is_set = .false.
    logical  :: unread = .false.
    logical  :: has_text = .false.
    character(len=:), allocatable :: text
  end type

  type(saved_value_t), allocatable :: saved(:)   ! Before a reading during the run (Save_Values)
  type(saved_value_t), allocatable :: start(:)   ! The reading of the set-up (Save_Start_Values)


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

  subroutine add_int(this,section,name,var,default,desc,allowed,required,during_run)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, intent(in) :: during_run   ! What a reading during the run does with the key: KEY_TAKEN, KEY_FIXED, ...
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
    this%params(n)%during_run = during_run

    this%params(n)%type_id = TYPE_INT
    this%params(n)%value%i => var

    read(default,*) var

  end subroutine


  !========================================================
  ! Add real parameter
  !========================================================

  subroutine add_real(this,section,name,var,default,desc,allowed,required,during_run)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, intent(in) :: during_run   ! What a reading during the run does with the key: KEY_TAKEN, KEY_FIXED, ...
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
    this%params(n)%during_run = during_run

    this%params(n)%type_id = TYPE_REAL
    this%params(n)%value%r => var

    read(default,*) var

  end subroutine add_real


  !========================================================
  ! Add logical parameter
  !========================================================

  subroutine add_logical(this,section,name,var,default,desc,allowed,required,during_run)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, intent(in) :: during_run   ! What a reading during the run does with the key: KEY_TAKEN, KEY_FIXED, ...
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
    this%params(n)%during_run = during_run

    this%params(n)%type_id = TYPE_LOG
    this%params(n)%value%l => var

    read(default,*) var

  end subroutine add_logical


  !========================================================
  ! Add string parameter
  !========================================================

  subroutine add_string(this,section,name,var,default,desc,allowed,required,during_run)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, intent(in) :: during_run   ! What a reading during the run does with the key: KEY_TAKEN, KEY_FIXED, ...
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
    this%params(n)%during_run = during_run

    this%params(n)%type_id = TYPE_STR
    this%params(n)%value%s => var

    var = default

  end subroutine add_string


  !========================================================
  ! Integer array
  !========================================================

  subroutine add_int_array(this,section,name,var,default,desc,allowed,required,during_run)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, intent(in) :: during_run   ! What a reading during the run does with the key: KEY_TAKEN, KEY_FIXED, ...
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
    this%params(n)%during_run = during_run

    this%params(n)%type_id = TYPE_INT
    this%params(n)%value%iarr => var

    read(default,*) defval
    var(:) = defval

  end subroutine


  !========================================================
  ! Real array
  !========================================================

  subroutine add_real_array(this,section,name,var,default,desc,allowed,required,during_run)

    implicit none

    class(registry_t), intent(inout) :: this
    character(*), intent(in) :: section,name,default,desc,allowed
    logical, intent(in) :: required
    integer, intent(in) :: during_run   ! What a reading during the run does with the key: KEY_TAKEN, KEY_FIXED, ...
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
    this%params(n)%during_run = during_run

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
  !> another form stops the run as well, so that no rule is skipped in silence, and so does a key that the registry holds
  !> with no class for a reading during the run (KEY_TAKEN, KEY_FIXED, KEY_FIXED_OWN, KEY_START). So does a real value that is not finite,
  !> with a rule or without, before the rule is checked: every comparison with a NaN is false, and an infinity passes a
  !> bound such as > 0 (gfortran reads a number too large for a real, such as 1e999, as an infinity). The test reads the
  !> bits (Not_Finite). Before every rule, a key of input.ini whose value does not read as the type of its parameter
  !> (Load_Ini marks it: a value that does not convert, NaN and Infinity included, an integer written with a period or
  !> an exponent, a list with more values than the parameter holds, no value) stops the run, and so does such a value of
  !> levels, which the reader takes outside the registry (Scan_Ini): without the stop, the run would go on with the
  !> default, or with the value of the previous reading, without a word.
  !> The message gives the file, the section, the key and the value as written.
  !> What an error does is up to the caller: at set-up it stops the run (Check_Input in Wrap_Setup); in a reading during
  !> the run it is an error of that reading, which ini-error decides (Read_Inifile_Runtime).
  function Validate_Registry() result(out)
    use MOSE_Finite_m, only: Not_Finite
    implicit none

    character(len=1024) :: out
    integer :: i, k
    real(R8) :: val
    character(len=16) :: nth

    out=""

    if (allocated(reg%load_error)) then
      out = reg%load_error
      return
    end if

    do i=1,reg%size

      if (reg%params(i)%during_run < KEY_TAKEN .or. reg%params(i)%during_run > KEY_START) then
        out = "[ERROR] "//trim(reg%params(i)%name)//": the registry gives it no class for a reading during the run"
        return
      end if

      if (reg%params(i)%unread) then
        out = Unread_Message(i)
        return
      end if

      if (reg%params(i)%required .and. .not. reg%params(i)%is_set) then
        out = "[ERROR] Required parameter not set: "//trim(reg%params(i)%name)
        return
      end if

      ! A real value, or a value of a list of reals, that is not finite
      if (associated(reg%params(i)%value%r)) then
        if (Not_Finite(reg%params(i)%value%r)) then
          out = "[ERROR] "//Key_And_Value(i)//": the value is not a finite real number"
          return
        end if
      else if (associated(reg%params(i)%value%rarr)) then
        do k = 1, size(reg%params(i)%value%rarr)
          if (Not_Finite(reg%params(i)%value%rarr(k))) then
            write(nth,'(I0)') k
            out = "[ERROR] "//Key_And_Value(i)//": value "//trim(nth)//" of the list is not a finite real number"
            return
          end if
        end do
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

        ! The reader takes a logical value only (a list-directed read: T or F, after an optional period, starts the value;
        ! any other value is a value that does not read, given above): nothing is left to check, if the rule says no more
        ! than that
        if (.not. logical_rule(reg%params(i)%allowed)) then
          out = "[ERROR] "//trim(reg%params(i)%name)//": the registry does not check its rule "//trim(reg%params(i)%allowed)
          return
        end if

      case default
        cycle

      end select

    end do

  end function Validate_Registry


  !> Keeps a copy of the value of every parameter of the registry, whatever its type, and of the fields that Load_Ini
  !> writes (is_set, unread, text): Restore_Values puts back the state of every parameter as it was at this call. Used
  !> when input.ini is read again during the run (Read_Inifile_Runtime), before Load_Ini.
  subroutine Save_Values()

  implicit none

  integer :: i

  if (allocated(saved)) then
    if (size(saved) /= reg%size) deallocate(saved)
  end if
  if (.not. allocated(saved)) allocate(saved(reg%size))
  do i = 1, reg%size
    call Copy_Value(i, saved(i))
  end do

  end subroutine Save_Values


  !> Puts back the value of every parameter of the registry, and the fields that Load_Ini writes, as Save_Values kept
  !> them: the variables of the parameters (the registry points at them) hold again what they held at that call
  subroutine Restore_Values()

  implicit none

  integer :: i

  if (.not. allocated(saved)) return
  do i = 1, min(reg%size, size(saved))
    call Put_Back(i)
  end do

  end subroutine Restore_Values


  !> Keeps a copy of the values that the reading of input.ini at set-up gives every parameter (Read_Inifile, right after
  !> Load_Ini, before the set-up makes its choices and writes some of them over: the name of a probe file, the shock
  !> detector that a Riemann solver carries): the values against which a reading during the run checks the keys fixed for
  !> the run (Fixed_Key_Changed)
  subroutine Save_Start_Values()

  implicit none

  integer :: i

  if (allocated(start)) deallocate(start)
  allocate(start(reg%size))
  do i = 1, reg%size
    call Copy_Value(i, start(i))
  end do

  end subroutine Save_Start_Values


  !> The error of the first key fixed for the run (KEY_FIXED) whose value in input.ini read again during the run is not
  !> the one that the reading of the set-up gave it (Save_Start_Values), or an empty string. Only the keys that the file
  !> holds count (Load_Ini keeps their text): a key that the file leaves out keeps its value, as every key does. The values
  !> are compared as the type of the key, so that 1e-5 and 1.0e-05 are the same value: integers and logicals for equality,
  !> reals with two ordered comparisons (Validate_Registry has found them finite), strings without their trailing blanks,
  !> lists value by value. The error gives the key and the value written, and the value of the set-up: as the file of the
  !> set-up wrote it, or, for a key that it left out, its default. Called by Read_Inifile_Runtime after Load_Ini and
  !> Validate_Registry; a key whose model checks it (KEY_FIXED_OWN) is left to the model
  function Fixed_Key_Changed() result(out)

  implicit none

  character(len=:), allocatable :: out
  character(len=:), allocatable :: name
  integer :: i

  out = ""
  if (.not. allocated(start)) return
  do i = 1, min(reg%size, size(start))
    if (reg%params(i)%during_run /= KEY_FIXED) cycle
    if (.not. allocated(reg%params(i)%text)) cycle
    if (.not. Value_Differs(i, start(i))) cycle
    name = trim(reg%params(i)%name)
    if (start(i)%has_text) then
      out = "[ERROR] "//Key_And_Value(i)//": the run has "//name//" = "//trim(adjustl(start(i)%text))//" from its set-up"
    else if (len_trim(reg%params(i)%default_str) > 0) then
      out = "[ERROR] "//Key_And_Value(i)//": the run started without "//name//" in input.ini (its default: "// &
            trim(reg%params(i)%default_str)//")"
    else
      out = "[ERROR] "//Key_And_Value(i)//": the run started without "//name//" in input.ini"
    end if
    out = out//", and "//name//" cannot change during the run"
    return
  end do

  end function Fixed_Key_Changed


  !> After a reading of input.ini during the run that the run takes (Read_Inifile_Runtime): the variables of the keys fixed
  !> for the run and of those of the start only hold again what they held before the reading (the copy of Save_Values), so
  !> that no code reads in them a value of the reading, and the values that the set-up wrote over stay
  subroutine Keep_Fixed_Values()

  implicit none

  integer :: i

  if (.not. allocated(saved)) return
  do i = 1, min(reg%size, size(saved))
    if (reg%params(i)%during_run /= KEY_TAKEN) call Put_Back(i)
  end do

  end subroutine Keep_Fixed_Values


  !> A copy of the value of parameter i, whatever its type, and of the fields that Load_Ini writes (is_set, unread, text)
  subroutine Copy_Value(i, copy)

  implicit none

  integer, intent(in) :: i
  type(saved_value_t), intent(inout) :: copy

  if (associated(reg%params(i)%value%i))    copy%i    = reg%params(i)%value%i
  if (associated(reg%params(i)%value%r))    copy%r    = reg%params(i)%value%r
  if (associated(reg%params(i)%value%l))    copy%l    = reg%params(i)%value%l
  if (associated(reg%params(i)%value%s))    copy%s    = reg%params(i)%value%s
  if (associated(reg%params(i)%value%iarr)) copy%iarr = reg%params(i)%value%iarr
  if (associated(reg%params(i)%value%rarr)) copy%rarr = reg%params(i)%value%rarr
  copy%is_set   = reg%params(i)%is_set
  copy%unread   = reg%params(i)%unread
  copy%has_text = allocated(reg%params(i)%text)
  if (copy%has_text) copy%text = reg%params(i)%text

  end subroutine Copy_Value


  !> Parameter i, its variable and the fields that Load_Ini writes, as Save_Values kept them
  subroutine Put_Back(i)

  implicit none

  integer, intent(in) :: i

  if (associated(reg%params(i)%value%i))    reg%params(i)%value%i    = saved(i)%i
  if (associated(reg%params(i)%value%r))    reg%params(i)%value%r    = saved(i)%r
  if (associated(reg%params(i)%value%l))    reg%params(i)%value%l    = saved(i)%l
  if (associated(reg%params(i)%value%s))    reg%params(i)%value%s    = saved(i)%s
  if (associated(reg%params(i)%value%iarr)) reg%params(i)%value%iarr = saved(i)%iarr
  if (associated(reg%params(i)%value%rarr)) reg%params(i)%value%rarr = saved(i)%rarr
  reg%params(i)%is_set = saved(i)%is_set
  reg%params(i)%unread = saved(i)%unread
  if (allocated(reg%params(i)%text)) deallocate(reg%params(i)%text)
  if (saved(i)%has_text) reg%params(i)%text = saved(i)%text

  end subroutine Put_Back


  !> Whether the value of parameter i differs from the copy, as the type of the parameter: integers and logicals for
  !> equality, reals with two ordered comparisons (an exact test, with no warning on real equality), strings without their
  !> trailing blanks, lists value by value
  logical function Value_Differs(i, copy)

  implicit none

  integer, intent(in) :: i
  type(saved_value_t), intent(in) :: copy

  if (associated(reg%params(i)%value%i)) then
    Value_Differs = reg%params(i)%value%i /= copy%i
  else if (associated(reg%params(i)%value%r)) then
    Value_Differs = reg%params(i)%value%r < copy%r .or. reg%params(i)%value%r > copy%r
  else if (associated(reg%params(i)%value%l)) then
    Value_Differs = reg%params(i)%value%l .neqv. copy%l
  else if (associated(reg%params(i)%value%s)) then
    Value_Differs = trim(reg%params(i)%value%s) /= trim(copy%s)
  else if (associated(reg%params(i)%value%iarr)) then
    Value_Differs = any(reg%params(i)%value%iarr /= copy%iarr)
  else if (associated(reg%params(i)%value%rarr)) then
    Value_Differs = any(reg%params(i)%value%rarr < copy%rarr .or. reg%params(i)%value%rarr > copy%rarr)
  else
    Value_Differs = .false.
  end if

  end function Value_Differs


  !> Whether a reading of input.ini taken by the run gave the key name of the section section a value of the type of its
  !> parameter (is_set, which Load_Ini sets and never clears): false for a key that every such reading left out, whose
  !> parameter keeps its default, and for a key not in the registry; set once, true also after a reading that leaves it out
  logical function Key_Is_Set(section, name)

  implicit none

  character(len=*), intent(in) :: section, name
  integer :: i

  Key_Is_Set = .false.
  do i = 1, reg%size
    if (reg%params(i)%section == section .and. reg%params(i)%name == name) then
      Key_Is_Set = reg%params(i)%is_set
      return
    end if
  end do

  end function Key_Is_Set


  !> The key i of input.ini, as the messages of Validate_Registry give it: the file, the section, the key and the value as
  !> written, when the file holds one
  function Key_And_Value(i) result(out)

  implicit none

  integer, intent(in) :: i
  character(len=:), allocatable :: out

  out = "input.ini, ["//trim(reg%params(i)%section)//"] "//trim(reg%params(i)%name)//" ="
  if (allocated(reg%params(i)%text)) out = out//" "//trim(reg%params(i)%text)

  end function Key_And_Value


  !> The error of a key of input.ini whose value does not read as the type of its parameter (Load_Ini): no value, or a
  !> value not of its type (a list: more values than the parameter holds, or a value not of its type)
  function Unread_Message(i) result(out)

  implicit none

  integer, intent(in) :: i
  character(len=:), allocatable :: out
  character(len=16) :: n

  out = "[ERROR] "//Key_And_Value(i)//": the value "
  if (.not. allocated(reg%params(i)%text)) then
    out = out//"is empty"
  else if (len_trim(reg%params(i)%text) == 0) then
    out = out//"is empty"
  else if (associated(reg%params(i)%value%iarr)) then
    write(n,'(I0)') size(reg%params(i)%value%iarr)
    out = out//"is not a list of at most "//trim(n)//" integers"
  else if (associated(reg%params(i)%value%rarr)) then
    write(n,'(I0)') size(reg%params(i)%value%rarr)
    out = out//"is not a list of at most "//trim(n)//" finite real numbers"
  else
    select case (reg%params(i)%type_id)
    case (TYPE_INT)
      out = out//"is not an integer"
    case (TYPE_REAL)
      out = out//"is not a finite real number"
    case (TYPE_LOG)
      out = out//"is not a logical (true or false)"
    case default
      out = out//"cannot be read"
    end select
  end if

  end function Unread_Message


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
    write(unit,'(A)') "MOSE reads `input.ini` again during the run, every `ini-diter` iterations. The column During the run "// &
      "says what such a reading does with each key: *taken*, the run takes the value written; *fixed*, the run keeps the "// &
      "value that it read at its start, and another value written is an error of the reading, which `ini-error` decides; "// &
      "*start only*, the key serves the start of the run, and a reading leaves it as it was."

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
        write(unit,'(A)') "| Parameter | Default | Allowed | Required | During the run | Description |"
        write(unit,'(A)') "|-----------|---------|---------|----------|----------------|-------------|"
      end if

      write(unit,'(A)') "| "//trim(this%params(i)%name)// &
          " | "//trim(this%params(i)%default_str)// &
          " | "//trim(this%params(i)%allowed)// &
          " | "//merge("yes"," no",this%params(i)%required)// &
          " | "//During_Run_Word(this%params(i)%during_run)// &
          " | "//trim(this%params(i)%description)//" |"

    end do

    close(unit)

  end subroutine generate_markdown


  !> The word of the column During the run of docs/user/registry.md for a class of key
  function During_Run_Word(during_run) result(word)

    implicit none

    integer, intent(in) :: during_run
    character(len=:), allocatable :: word

    select case (during_run)
    case (KEY_TAKEN)
      word = "taken"
    case (KEY_FIXED, KEY_FIXED_OWN)
      word = "fixed"
    case (KEY_START)
      word = "start only"
    case default
      word = ""
    end select

  end function During_Run_Word


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