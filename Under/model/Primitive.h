#pragma once

#define LEFT_SELECTED        0x000001
#define RIGHT_SELECTED       0x000002
#define SELECTED             0x000003

#define DELETED              0x000004
#define HIDDEN_TAG             0x000008


class CPrimitive
{
public:
	CPrimitive();
	~CPrimitive();

	void SetLeftSelected(bool bSelect = true)
	{
		SetTag(LEFT_SELECTED, bSelect);
	}

	bool IsLeftSelected() const
	{
		return IsTagSet(LEFT_SELECTED);
	}

	void SetRightSelected(bool bSelect = true)
	{
		SetTag(RIGHT_SELECTED, bSelect);
	}

	bool IsRightSelected() const
	{
		return IsTagSet(RIGHT_SELECTED);
	}

	void SetSelected(bool bSelect = true)
	{
		SetTag(SELECTED, bSelect);
	}

	bool IsSelected() const
	{
		return IsTagSet(SELECTED);
	}

	void SetDeleted(bool bDeleted = true)
	{
		SetTag(DELETED, bDeleted);
	}

	bool IsDeleted() const
	{
		return IsTagSet(DELETED);
	}

	void SetShown(bool show = true)
	{
		SetTag(HIDDEN_TAG, !show);
	}

	bool IsShown() const
	{
		return !IsTagSet(HIDDEN_TAG);
	}

private:
	void SetTag(unsigned int tag, bool set)
	{
		if (set)
		{
			m_tag |= tag;
		}
		else
		{
			m_tag &= (~tag);
		}
	}

	bool IsTagSet(unsigned int tag) const
	{
		return (m_tag & tag) != 0;
	}

private:
	unsigned int m_tag;
};

