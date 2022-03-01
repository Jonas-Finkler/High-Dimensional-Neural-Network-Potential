! Copyright (C) 2020 Jonas A. Finkler
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


module argumentParser
    use precision

    type argument
        character(len=1) :: name
        logical :: required
        logical :: wasFound
    end type argument

    type, extends(argument) :: stringArgument
        character(:), allocatable :: default
        character(:), allocatable :: value
    end type stringArgument

    type, extends(argument) :: intArgument
        integer :: default
        integer :: value
    end type intArgument

    type, extends(argument) :: realArgument
        real(dp) :: default
        real(dp) :: value
    end type realArgument

    type, extends(argument) :: logicalArgument
        logical :: default
        logical :: value
    end type

    type argumentListElement
        class(argument), pointer :: el
    end type argumentListElement

    type argumentList
        character(:), allocatable :: programName
        character(:), allocatable :: helpText
        character(:), allocatable :: version
        integer :: size
        integer :: nArgs
        type(argumentListElement), allocatable :: list(:)
    end type argumentList

contains

    subroutine initArgumentParser(arglist, programName, helpText, version)
        type(argumentList), intent(out) :: argList
        character(len=*) :: programName
        character(len=*) :: helpText
        character(len=*) :: version

        argList%programName = programName
        argList%helpText = helpText
        argList%version = version

        argList%nArgs = 0
        argList%size = 1
        allocate(argList%list(argList%size))

    end subroutine initArgumentParser

    subroutine addArgument(argList, arg)
        type(argumentList), intent(inout) :: argList
        class(argument), target, intent(inout) :: arg
        type(argumentListElement), allocatable :: tmpList(:)

        select type(arg)
            type is (argument)
                print*, 'ERROR: Only subtypes of argument can be passed to addArgument'
                stop
        end select

        arg%wasFound = .false.

        if (argList%nArgs >= argList%size) then
            allocate(tmpList(argList%nArgs))
            tmpList = argList%list
            deallocate(argList%list)
            allocate(argList%list(argList%size*2))
            argList%list(:argList%size) = tmpList(:)
            argList%size = argList%size * 2
        end if
        argList%nArgs = argList%nArgs + 1
        argList%list(argList%nArgs)%el => arg

    end subroutine addArgument

    subroutine stopProgram(argList)
        type(argumentList), intent(in) :: argList
        print*, 'Use "' // argList%programName // ' -h" to find help.'
        stop
    end subroutine

    subroutine parseArguments(argList)
        type(argumentList), intent(inout) :: argList
        integer :: nargs
        integer :: iarg
        character(len=200) arg
        integer :: arglen, argstatus
        integer :: i
        class(argument), pointer :: larg
        character(len=200) :: tmpstring
        logical :: argExists

        nargs = command_argument_count()

        ! todo: Maybe allow spaced between argument and value
        do iarg=1,nargs
            call get_command_argument(iarg, arg, arglen, argstatus)
            if (arg(1:1) /= '-') then
                ! todo: print help message here
                print*, "Arguments must beginn with '-'"
                call stopProgram(argList)
                stop
            end if
            if (arg(1:2) == '-h' .or. arg(1:3) == '--h') then
                print*, argList%helpText
                stop
            end if
            if (arg(1:2) == '-v' .or. arg(1:3) == '--v') then
                print*, 'Version: ' // argList%version
                stop
            end if
            argExists = .false.
            do i=1, argList%nArgs
                larg => argList%list(i)%el
                if (arg(2:2) == larg%name) then
                    if (larg%wasFound) then
                        print*, 'Error: Argument ', larg%name, ' was specified twice'
                        call stopProgram(argList)
                    end if
                    larg%wasFound = .true.
                    select type(larg)
                        type is (stringArgument)
                            read(arg(3:),*) tmpstring
                            larg%value = trim(tmpstring)
                            ! print*, 'Found Argumet ', larg%name, ' =', larg%value
                        type is (intArgument)
                            read(arg(3:),*) larg%value
                            ! print*, 'Found Argumet ', larg%name, ' =', larg%value
                        type is (realArgument)
                            read(arg(3:),*) larg%value
                            ! print*, 'Found Argumet ', larg%name, ' =', larg%value
                        type is (logicalArgument)
                            read(arg(3:),*) larg%value
                            ! print*, 'Found Argumet ', larg%name, ' =', larg%value
                    end select
                    argExists = .true.
                    exit
                end if

            end do
            if (.not. argExists) then
                print*, 'Error: Unknown argument ', arg
                call stopProgram(argList)
            end if
        end do

        do i=1, argList%nArgs
            larg => argList%list(i)%el
            if (.not. larg%wasFound) then
                if (larg%required) then
                    print*, 'ERROR: Argument -', larg%name, ' is required'
                    call stopProgram(argList)
                else
                    select type(larg)
                        type is (stringArgument)
                            larg%value = larg%default
                        type is (intArgument)
                            larg%value = larg%default
                        type is (realArgument)
                            larg%value = larg%default
                        type is (logicalArgument)
                            larg%value = larg%default
                    end select
                end if
            end if
        end do

    end subroutine parseArguments


end module argumentParser