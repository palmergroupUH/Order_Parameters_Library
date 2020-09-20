module linked_list
    use system
    implicit none 

    private

    type linkedlist_type
        integer :: element
        type(linkedlist_type), pointer :: next
    end type

    public :: linkedlist_type, &
            & LinkedlistAdd, &
            & LinkedListPrint, & 
            & LinkedListCount, &
            & Linkedlist2Array 

contains

    function LinkedlistAdd(head,newelement)
        implicit none
        integer,intent(in) :: newelement
        type(linkedlist_type),pointer,intent(in) :: head
        type(linkedlist_type),pointer :: new,LinkedlistAdd

        allocate(new)

        new%element = newelement

        new%next => head

        LinkedlistAdd => new

        end function

    subroutine LinkedListPrint(head)
        implicit none
        type( linkedlist_type ), pointer,intent(in) :: head
        type( linkedlist_type ), pointer :: temp

        temp => head

        do while ( associated(temp) )

            print*, temp%element

            temp => temp%next

        end do

        end subroutine

    integer function LinkedListCount(head)
        implicit none
        type(linkedlist_type), pointer,intent(in) :: head
        type(linkedlist_type), pointer :: temp

        temp => head

        LinkedListCount = 0

        do while ( associated(temp))

            LinkedListCount = LinkedListCount + 1

            temp => temp%next

        end do

        end function

    function Linkedlist2Array(head,total_num) result(array)
        implicit none
        integer,intent(in) :: total_num
        type(linkedlist_type), pointer,intent(in) :: head
        type(linkedlist_type), pointer :: temp
        integer :: counter
        integer,dimension(1:total_num) :: array

        temp => head

        counter = 0

        do while ( associated(temp))

            counter = counter + 1

            array(counter) = temp%element

            temp => temp%next

        end do

        end function

end module 
